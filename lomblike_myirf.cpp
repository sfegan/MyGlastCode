//-*-mode:c++; mode:font-lock;-*-

// lomelike.cpp - Stephen Fegan 
//              - sfegan@llr.in2p3.fr
//              - 2010-07-07
//
// Apply "Lomb"-like Likelihood test for oscillatory model
//
// $Id: lomblike_myirf.cpp 2819 2011-05-19 15:19:43Z sfegan $

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"
#include "VERITAS/VSFileUtility.hpp"
#include "VERITAS/VSDataConverter.hpp"

#include "FT1ROIFilter.hpp"
#include "FT2Exp.hpp"
#include "IRF.hpp"
#include "PLIntegratedEA.hpp"
#include "GTI.hpp"
#include "FITS.hpp"
#include "Util.hpp"
#include "MyMinuit.hpp"
#include "Accumulator.hpp"
#include "MyMinuit.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

// ****************************************************************************
//
// PERIODGRAM
//
// ****************************************************************************

class OscillatoryModel: public Optimizable
{
public:
  class Wave
  {
  public:
    Wave(): 
      freq(0), amp(0), phase(0), amp_free(false), phase_free(false) { }
    Wave(double f, double a=0, double p=0, bool afree=true, bool pfree=true):
      freq(f), amp(a), phase(p), amp_free(afree), phase_free(pfree) { }
    double freq;
    double amp;
    double phase;
    bool amp_free;
    bool phase_free;
  };

  OscillatoryModel(const std::vector<Wave>& waves,
		   const FT2Exp& ft2exp, const std::vector<FT1AEff>& events);
  virtual ~OscillatoryModel();

  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);

  double aNorm() const { return m_anorm; }
  const std::vector<Wave>& waves() { return m_wave; }

private:
  void setP(const std::vector<double>& p)
  {
    unsigned nfree = 0;
    for(unsigned iwave=0;iwave<m_wave.size();iwave++)
      {
	if(m_wave[iwave].amp_free)
	  {
	    if(m_wave[iwave].amp!=p[nfree])
	      m_wave[iwave].amp=p[nfree];
	    nfree++;
	  }
	if(m_wave[iwave].phase_free)
	  {
	    if(m_wave[iwave].phase!=p[nfree])
	      m_wave[iwave].phase=p[nfree],
		m_cphi[iwave]=std::cos(m_wave[iwave].phase),
		m_sphi[iwave]=std::sin(m_wave[iwave].phase);
	    nfree++;
	  }
      }
  }

  inline double evPrFactor(unsigned iwave, unsigned ievent) const;
  inline double dEvPrFactor_dPhi(unsigned iwave, unsigned ievent) const;
  inline double evPr(unsigned ievent) const;

  inline double nPredFactor(unsigned iwave) const;
  inline double dNPredFactor_dPhi(unsigned iwave) const;

  double dfdAmp(unsigned iwave) const;
  double dfdPhi(unsigned iwave) const;
  
  std::vector<Wave> m_wave;
  const FT2Exp& m_ft2exp;
  const std::vector<FT1AEff>& m_events;
  double m_exposure;
  double m_anorm;

  std::vector<double> m_cphi;
  std::vector<double> m_sphi;

  std::vector<double> m_int_ecwt;
  std::vector<double> m_int_eswt;
  std::vector<std::vector<double> > m_cwti;
  std::vector<std::vector<double> > m_swti;
};

OscillatoryModel::
OscillatoryModel(const std::vector<Wave>& waves,
		 const FT2Exp& ft2exp, const std::vector<FT1AEff>& events):
  Optimizable(), m_wave(waves), m_ft2exp(ft2exp), m_events(events),
  m_exposure(), m_anorm(),
  m_cphi(m_wave.size()), 
  m_sphi(m_wave.size()),
  m_int_ecwt(m_wave.size()),
  m_int_eswt(m_wave.size()),
  m_cwti(m_wave.size(),std::vector<double>(m_events.size())),
  m_swti(m_wave.size(),std::vector<double>(m_events.size()))
{
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      m_wave[iwave].freq
	=fabs(m_wave[iwave].freq);
      m_wave[iwave].phase
	=fmod(fmod(m_wave[iwave].phase,2.0*M_PI)+2.0*M_PI,2.0*M_PI);
      if(m_wave[iwave].freq==0.0)
	{
	  m_wave[iwave].phase=0;
	  m_wave[iwave].phase_free=false;
	}
    }
  
  for(unsigned ie=0;ie<m_ft2exp.nEInt();ie++)
    m_exposure += m_ft2exp.eInt(ie).back().stop; 
  m_anorm = double(m_events.size())/m_exposure;

  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      const double wi = M_PI*2.0*m_wave[iwave].freq;
      m_wave[iwave].freq = wi;

      m_cphi[iwave] = std::cos(m_wave[iwave].phase);
      m_sphi[iwave] = std::sin(m_wave[iwave].phase);

      if(wi == 0)
	{
	  m_int_ecwt[iwave] = m_exposure;
	  m_int_eswt[iwave] = 0;
	}
      else
	{
	  Accumulator ac;
	  Accumulator as;
	  for(unsigned jint=0;jint<m_ft2exp.nInt();jint++)
	    {
	      const double witj = wi*m_ft2exp.tInt()[jint].start;
	      double ej = 0;
	      for(unsigned ie=0;ie<m_ft2exp.nEInt();ie++)
		if(!m_ft2exp.eInt(ie).empty())
		  ej += m_ft2exp.eInt(ie)[jint].x2;
	      ac.add(ej*std::cos(witj));
	      as.add(ej*std::sin(witj));
	    }
	  m_int_ecwt[iwave] = ac.sum();
	  m_int_eswt[iwave] = as.sum();
	}
      
      if(wi == 0.0)
	for(unsigned jevent=0;jevent<m_events.size();jevent++)
	  {
	    m_cwti[iwave][jevent] = 1.0;
	    m_swti[iwave][jevent] = 0.0;
	  }
      else
	for(unsigned jevent=0;jevent<m_events.size();jevent++)
	  {
	    const double witj = wi*m_events[jevent].time;
	    m_cwti[iwave][jevent] = std::cos(witj);
	    m_swti[iwave][jevent] = std::sin(witj);
	  }
    }
}

OscillatoryModel::~OscillatoryModel()
{
  // nothing to see here
}

unsigned OscillatoryModel::numParam()
{
  unsigned nfree = 0;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      if(m_wave[iwave].amp_free)nfree++;
      if(m_wave[iwave].phase_free)nfree++;
    }
  return nfree;
}

std::string OscillatoryModel::pName(unsigned iparam)
{
  unsigned nfree = 0;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      if(m_wave[iwave].amp_free)nfree++;
      if(nfree>iparam)
	return std::string("AMP")+VSDataConverter::toString(iwave);
      if(m_wave[iwave].phase_free)nfree++;
      if(nfree>iparam)
	return std::string("PHI")+VSDataConverter::toString(iwave);
    }
  assert(0);
}

bool OscillatoryModel::canCalcDFDP()
{
  return true; //false;
}

inline double OscillatoryModel::
evPrFactor(unsigned iwave, unsigned ievent) const
{
  if(m_wave[iwave].freq==0.0)return 1;  
  else return (1.0 
	       + m_swti[iwave][ievent]*m_cphi[iwave]
	       - m_cwti[iwave][ievent]*m_sphi[iwave]);
}

inline double OscillatoryModel::
dEvPrFactor_dPhi(unsigned iwave, unsigned ievent) const
{
  if(m_wave[iwave].freq==0.0)return 0;  
  else return (- m_swti[iwave][ievent]*m_sphi[iwave] 
	       - m_cwti[iwave][ievent]*m_cphi[iwave]);
}

inline double OscillatoryModel::evPr(unsigned ievent) const
{
  Accumulator a;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    a.add(m_wave[iwave].amp*evPrFactor(iwave,ievent));
  return a.sum();
}

inline double OscillatoryModel::nPredFactor(unsigned iwave) const
{
  // If statement not strictly necessary - 2nd formula would do for both
  if(m_wave[iwave].freq==0.0)return m_exposure;
  else return (m_exposure 
	       + m_int_eswt[iwave]*m_cphi[iwave]
	       - m_int_ecwt[iwave]*m_sphi[iwave]);
}

inline double OscillatoryModel::dNPredFactor_dPhi(unsigned iwave) const
{
  // If statement not strictly necessary - 2nd formula would do for both
  if(m_wave[iwave].freq==0.0)return 0;
  else return (- m_int_eswt[iwave]*m_sphi[iwave]
	       - m_int_ecwt[iwave]*m_cphi[iwave]);
}


double OscillatoryModel::f(const std::vector<double>& p)
{
  setP(p);
  Accumulator a;

  // Data sum
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      double p = evPr(ievent);
      assert(p>0);
      if(p>0)a.add(std::log(evPr(ievent)));
    }
  
  // Model sum
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    a.add(-m_wave[iwave].amp*m_anorm*nPredFactor(iwave));

  return -a.sum();
}

double OscillatoryModel::dfdAmp(unsigned iwave) const
{
  Accumulator a;
 
  // Data sum
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      double p = evPr(ievent);
      if(p>0)a.add(evPrFactor(iwave,ievent)/p);
    }

  // Model sum
  a.add(-m_anorm*nPredFactor(iwave));

  return -a.sum();  
}

double OscillatoryModel::dfdPhi(unsigned iwave) const
{
  if(m_wave[iwave].freq==0.0)return 0;
 
  Accumulator a;
  // Data sum
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      double p = evPr(ievent);
      if(p>0)a.add(dEvPrFactor_dPhi(iwave,ievent)/p);
    }

  // Model sum
  a.add(-m_anorm*dNPredFactor_dPhi(iwave));

  return -a.sum()*m_wave[iwave].amp;
}

double OscillatoryModel::dfdp(const std::vector<double>& p, unsigned iparam)
{
  setP(p);
  unsigned nfree = 0;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      if(m_wave[iwave].amp_free)nfree++;
      if(nfree>iparam)return dfdAmp(iwave);
      if(m_wave[iwave].phase_free)nfree++;
      if(nfree>iparam)return dfdPhi(iwave);
    }
  assert(0);
}

double OscillatoryModel::p0(unsigned iparam)
{
  unsigned nfree = 0;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      if(m_wave[iwave].amp_free)nfree++;
      if(nfree>iparam)return m_wave[iwave].amp;
      if(m_wave[iwave].phase_free)nfree++;
      if(nfree>iparam)return m_wave[iwave].phase;
    }
  assert(0);
}

double OscillatoryModel::plo(unsigned iparam)
{
  unsigned nfree = 0;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      if(m_wave[iwave].amp_free)nfree++;
      if(nfree>iparam)return 0.0;
      if(m_wave[iwave].phase_free)nfree++;
      if(nfree>iparam)return -2.0*M_PI;
    }
  assert(0);
}

double OscillatoryModel::phi(unsigned iparam)
{
  unsigned nfree = 0;
  for(unsigned iwave=0;iwave<m_wave.size();iwave++)
    {
      if(m_wave[iwave].amp_free)nfree++;
      if(nfree>iparam)return 2.0;
      if(m_wave[iwave].phase_free)nfree++;
      if(nfree>iparam)return 4.0*M_PI;
    }
  assert(0);
}

class Scanner
{
public:
  class Scan
  {

  };

  Scanner(double scan0, double scanN, unsigned nscan, 
	  const FT2Exp& ft2exp, unsigned verbose, bool log_scan = false,
	  bool phase_scan=false, double phase_scan_freq=1.0):
    m_scan0(scan0), m_scanN(scanN), m_nscan(nscan), m_ft2exp(ft2exp),
    m_verbose(verbose), m_log_scan(log_scan),
    m_phase_scan(phase_scan), m_phase_scan_freq(phase_scan_freq),
    m_null_logl(), m_null_amp(), m_null_damp() { }
    
  void scan(const std::vector<FT1AEff>& events);

private:
  double m_scan0;
  double m_scanN;
  unsigned m_nscan;
  const FT2Exp& m_ft2exp;
  unsigned m_verbose;

  bool m_log_scan;
  bool m_phase_scan;
  double m_phase_scan_freq;

  double m_null_logl;
  double m_null_amp;
  double m_null_damp;
};

void Scanner::scan(const std::vector<FT1AEff>& events)
{
  double exposure(0);
  for(unsigned ie=0;ie<m_ft2exp.nEInt();ie++)
    exposure += m_ft2exp.eInt(ie).back().stop;

  std::vector<OscillatoryModel::Wave> wave;
  wave.push_back(OscillatoryModel::Wave(0,0.9));

  if(1)
    {
      OscillatoryModel m(wave, m_ft2exp, events);
      MyMinuit minuit(&m,std::max(int(m_verbose)-4,-1));
      minuit.minimize(/*false,1e-3,0.5,"",true*/);
      m_null_logl = minuit.fVal();
      m_null_amp  = minuit.pVal(0) * m.aNorm();
      m_null_damp = minuit.pErr(0) * m.aNorm();
    }

  if(m_verbose>=2)
    {
      std::cout 
	<< "FIT: Const model A=" << m_null_amp << "+/-" << m_null_damp << ' '
	<< "LogL=" << m_null_logl << '\n';
    }

  wave.push_back(OscillatoryModel::Wave(m_phase_scan_freq,0.1,0,
					true,!m_phase_scan));

  for(unsigned iscan=0;iscan<m_nscan;iscan++)
    {
      double xscan = m_scan0;
      if(m_nscan>1)
	{
	  double x0 = m_scan0;
	  double xN = m_scanN;
	  if(m_log_scan)x0 = log(x0), xN=log(xN);
	  xscan = x0+(xN-x0)/double(m_nscan-1)*double(iscan);
	  if(m_log_scan)xscan = exp(xscan);
	}

      if(m_phase_scan)wave[1].phase = xscan;
      else wave[1].freq = xscan;

      OscillatoryModel m(wave, m_ft2exp, events);
      MyMinuit minuit(&m,std::max(int(m_verbose)-4,-1));

      try
	{
	  minuit.minimize(/*false,1e-3,0.5,"",true*/);
	}
      catch(const std::string& x)
	{
	  std::cout << x;
	}
      
      if(m_verbose>=1)
	{
	  if(m_phase_scan)
	    std::cout << xscan/(2.0*M_PI) << ' ' 
		      << xscan/(2.0*M_PI)/m_phase_scan_freq/86400.0 << ' ';
	  else
	    std::cout << xscan*86400.0 << ' ' 
		      << 1.0/(86400.0*xscan) << ' ';
	  std::cout 
	    << minuit.pVal(1)*m.aNorm() << ' ' 
	    << minuit.pErr(1)*m.aNorm() << ' '
	    << m_null_logl-minuit.fVal() << ' ';
	  if(!m_phase_scan)
	    std::cout 
	      << minuit.pVal(2)/(2.0*M_PI) << ' ' 
	      << minuit.pErr(2)/(2.0*M_PI) << ' ';
	  std::cout 
	    << minuit.pVal(0)*m.aNorm() << ' ' 
	    << minuit.pErr(0)*m.aNorm() << '\n';
	}
    }
}

// ****************************************************************************
//
// MAIN
//
// ****************************************************************************

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] ft1_file ft2_file irfs"
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);

  std::string command_line(*argv);
  for(int iarg=1;iarg<argc;iarg++)
    {
      command_line+=" ";
      command_line+=argv[iarg];
    }

  // --------------------------------------------------------------------------
  // PROCESS OPTIONS
  // --------------------------------------------------------------------------

  VSOptions options(argc, argv, true);
  
  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  int verbose = 1;
  if(options.find("q","Print no messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 0;
  if(options.find("v","Print really verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 2;
  if(options.find("vv","Print extremely verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 3;
  if(options.find("vvv","Print extremely verbose messages during analysis "
		  "and enable verbose output from the minimizer.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 4;

  // ---------------------------- SCANNING OPTIONS ----------------------------

  bool phase_scan = false;
  if(options.find("scan_phase","Do scan in phase (at fixed frequency) rather "
		  "than normal frequency/phase scan.","scan")
     !=VSOptions::FS_NOT_FOUND)
    phase_scan = 2;

  double scan0 = 1.0/100.0;
  double scanN = 1.0/1.0;
  unsigned nscan = 100;;
  double phase_scan_freq = 1.0/1.0;

  if(phase_scan)scan0=0, scanN=1.0;

  options.findWithValue("scan0", scan0,
			"Lowest frequency in scan [1/day] or lowest phase "
			"if doing a phase scan [cycle].","scan");
  options.findWithValue("scanN", scanN,
			"Highest frequency in scan [1/day] or highest phase "
			"if doing a phase scan [cycle].","scan");
  options.findWithValue("nscan", nscan,
			"Number of frequencies (or periods) in scan.","scan");

  bool log_scan = false;
  if(options.find("scan_log","Scan frequencies (or phases) logarithmically.",
		  "scan")
     !=VSOptions::FS_NOT_FOUND)
    log_scan = true;

  options.findWithValue("phase_scan_freq", phase_scan_freq,
			"Fixed frequency for phase scan [1/day].", "scan");

  if(phase_scan)scan0*=2.0*M_PI, scanN*=2.0*M_PI;
  else scan0 /= 86400.0, scanN /= 86400.0;
  phase_scan_freq /= 86400.0;

  // ---------------------------- EVENT CUT OPTIONS ---------------------------

  bool has_energy_cuts = false;
  std::pair<double,double> energy_cuts;
  bool has_position_cuts = false;
  triple<double,double,double> position_cuts;
  if(options.findWithValue("erange", energy_cuts,
			   "Energy range to consider in analysis (specifies "
			   "cuts on energy in FT1 file). If not provided, "
			   "the full range of energies in the FT1 file is "
			   "used. Specify as emin/emax in MeV.",
			   "cuts") != VSOptions::FS_NOT_FOUND)
    has_energy_cuts = true;
  
  if(options.findWithValue("roi", position_cuts,
			   "Circular region to consider in analysis "
			   "(specifies cuts on event position in FT1 file). "
			   "If not provided, the full ROI is used. Specify "
			   "as ra/dec/radus in degrees.",
			   "cuts") != VSOptions::FS_NOT_FOUND)
    {
      has_position_cuts = true;
      position_cuts.first  = d2r(position_cuts.first);
      position_cuts.second = d2r(position_cuts.second);
      position_cuts.third  = d2r(position_cuts.third);
    }

  // -------------------------------- EXPOSURE --------------------------------

  double gamma = 2.0;
  options.findWithValue("gamma", gamma,
			"Spectral index to use in weighting of effective "
			"area.", "exposure");

  bool no_livetime_correction = false;
  if(options.find("no_lt_correction","Do not apply livetime correction (if it "
		  "available in the IRFs)",
		  "exposure")
     !=VSOptions::FS_NOT_FOUND)
    no_livetime_correction = true;

  std::string day_exp_lc_fn = std::string();
  options.findWithValue("day_exp_lc", day_exp_lc_fn,
			"Specify the name of a file in which to write the "
			"lightcurve of daily exposure values. This is mostly "
			"meant as a diagnostics that can be compared with a "
			"curve produced by gtexposure.","exposure");

  std::string caldb = "$CALDB";
  options.findWithValue("caldb", caldb,
			"Base directory of calibration database.","exposure");

  // --------------------------------------------------------------------------
  // FINISH OPTIONS PROCESSING
  // --------------------------------------------------------------------------

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
	if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  int arg_req = 3;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string ft1_fn(*argv);
  argc--, argv++;

  std::string ft2_fn(*argv);
  argc--, argv++;

  std::string irf(*argv);
  argc--, argv++;  

  // --------------------------------------------------------------------------
  // GET ROI CUTS
  // --------------------------------------------------------------------------

  FITSHeader ft1head;
  if(verbose>=3)
    std::cout << "Load: FT1 header - " << ft1_fn << '\n';
  ft1head.loadFromFITS(ft1_fn,FT1::tableName());
  
  if(!ft1head.has("NDSKEYS"))
    {
      std::cerr << "Fatal: " << ft1_fn << " does not have ROI cuts\n";
      std::exit(EXIT_FAILURE);
    }

  double roi_emin(0);
  double roi_emax(0);
  double roi_ra(0);
  double roi_dec(0);
  double roi_roi(0);

  unsigned nds(0);
  ft1head.getAs(nds,"NDSKEYS");
  for(unsigned ids=0;ids<nds;ids++)
    {
      std::string type = ft1head.get(ft1head.nk("DSTYP",ids+1));
      if(type == "ENERGY")
	{
	  std::string val = ft1head.get(ft1head.nk("DSVAL",ids+1));
	  std::string::size_type ic = val.find(':');
	  if(ic != val.npos)
	    {
	      VSDataConverter::fromString(roi_emin,val.substr(0,ic));
	      VSDataConverter::fromString(roi_emax,val.substr(ic+1));
	    }
	}
      else if(type == "POS(RA,DEC)")
	{
	  std::string val = ft1head.get(ft1head.nk("DSVAL",ids+1));
	  if(val.substr(0,7)=="CIRCLE(")
	    {
	      std::vector<double> dv;
	      VSDataConverter::fromString(dv,val.substr(7,val.length()-8));
	      roi_ra  = d2r(dv.at(0));
	      roi_dec = d2r(dv.at(1));
	      roi_roi = d2r(dv.at(2));
	    }
	}
    }

  if(verbose>=3)
    {
      std::cout << "FT1: " << roi_emin << " <= E <= " << roi_emax << '\n';
      std::cout << "FT1: |(RA,Dec)-(" 
		<< r2d(roi_ra) << ',' << r2d(roi_dec)
		<< ")| <= " << r2d(roi_roi) << '\n';
    }

  double emin(roi_emin);
  double emax(roi_emax);
  if(has_energy_cuts)
    {
      if(energy_cuts.first<emin)
	{
	  std::cerr << "Minimum energy cut less than FT1 value: " 
		    << energy_cuts.first << " < " << roi_emin << '\n';
	  std::exit(EXIT_FAILURE);
	}

      if(energy_cuts.second>emax)
	{
	  std::cerr << "Maximum energy cut greater than FT1 value: " 
		    << energy_cuts.second << " > " << roi_emax << '\n';
	  std::exit(EXIT_FAILURE);
	}

      emin = energy_cuts.first;
      emax = energy_cuts.second;
    }
  double logemin(std::log10(emin));
  double logemax(std::log10(emax));

  double ra(roi_ra);
  double dec(roi_dec);
  double roi(roi_roi);
  if(has_position_cuts)
    {
      double d = position_cuts.third +
	sphere_dist(ra, dec, position_cuts.first, position_cuts.second);

      if(d > roi)
	{
	  std::cerr << "Desired ROI is outside FT1 file region: "
		    << r2d(d) << " > " << r2d(roi) << '\n';
	  std::exit(EXIT_FAILURE);
	}
      
      ra  = position_cuts.first;
      dec = position_cuts.second;
      roi = position_cuts.third;
    }

  // --------------------------------------------------------------------------
  // READ IN IRFS
  // --------------------------------------------------------------------------

  std::vector<IRFs*> irfs(2);
  for(unsigned iirf=0;iirf<irfs.size();iirf++)
    {
      irfs[iirf] = new IRFs;
      irfs[iirf]->loadAllIRFs(irf, iirf==0?true:false, !no_livetime_correction,
			      verbose>=1, caldb);
    }

  double eaint_ctheta0 = 0.2;
  double eaint_dctheta = 0.0025;

  std::vector<PLIntegratedEA*> eaint(irfs.size());
  for(unsigned iirf=0;iirf<irfs.size();iirf++)
    {
      if(irfs[iirf]->ea())eaint[iirf] = 
        new PLIntegratedEA(irfs[iirf], eaint_dctheta, eaint_ctheta0,
			   logemin, logemax, gamma);
    }

  // --------------------------------------------------------------------------
  // READ IN FT1 GTI
  // --------------------------------------------------------------------------

  GTIRange gti;
  if(verbose>=3)
    std::cout << "Load: GTI - " << ft1_fn << '\n';
  gti.loadGTIsFromFITS(ft1_fn);
  
  if(verbose>=2)
    {
      std::cout << "GTI: Total time " << lrint(gti.totalTime()) 
		<< " seconds in " << gti.nGTIs() << " periods\n";
    }

#if 0
  for(GTIRange::const_iterator igti=gti.begin();igti!=gti.end();igti++)
    std::cout << lrint(igti->t_start) << ' ' << lrint(igti->t_stop) << '\n';
#endif

  // --------------------------------------------------------------------------
  // READ IN FT2
  // --------------------------------------------------------------------------
  
  FT2Exp ft2exp(eaint, gti, ra, dec);
  FITSVectorDispatcher<FT2> ft2_dispatcher(&ft2exp);

  if(verbose>=2)
    std::cout << "Load: FT2 - " << ft2_fn << '\n';
  unsigned ft2_ntotal = 
    ft2_dispatcher.dispatchVector(ft2_fn, FT2::tableName());

  if(verbose>=2)
    {
      std::cout
	<< "FT2: " << ft2_ntotal << " entries, of which " 
	<< ft2exp.nInt() << " overlap with GTI\n";
      std::cout
	<< "FT2: Total time exposure: " << lrint(ft2exp.tInt().back().x2)
	<< " seconds\n";
      for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
	std::cout
	  << "FT2: Total exposure (EC=" << iea << "): " 
	  << ft2exp.eInt(iea).back().stop << " cm^2 seconds\n";
    }

  if(!day_exp_lc_fn.empty())
    {
      const std::vector<std::pair<unsigned,double> >& lc(ft2exp.getDayExpLC());
      if(verbose>=2)
	std::cout
	  << "FT2: writing exposure lightcurve to: " << day_exp_lc_fn << '\n';
      std::ofstream str(day_exp_lc_fn.c_str());
      for(unsigned ilc=0;ilc<lc.size();ilc++)
	str << lc[ilc].first << ' ' << lc[ilc].second << '\n';
    }

  // --------------------------------------------------------------------------
  // READ IN FT1
  // --------------------------------------------------------------------------

  FT1ROIFilter ft1(ft2exp);
  if(has_position_cuts)ft1.setROICuts(ra, dec, roi);
  if(has_energy_cuts)ft1.setEnergyCuts(emin, emax);
  FITSVectorDispatcher<FT1> ft1_dispatcher(&ft1);
  
  if(verbose>=2)
    std::cout << "Load: FT1 - " << ft1_fn << '\n';
  unsigned ft1_ntotal = 
    ft1_dispatcher.dispatchVector(ft1_fn, FT1::tableName());
  
  if(verbose>=2)
    {
      std:: cout
	<< "FT1: Found " << ft1_ntotal << " events, kept " 
	<< ft1.events.size() << " (";
      for(unsigned ict=0;ict<ft1.nCT();ict++)
	{
	  if(ict!=0)std::cout << ", ";
	  std::cout << "CT=" << ict << ":" << ft1.nEvents(ict);
	}
      std::cout	<< ")\n";
      if(ft1.nCutNoFT2())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutNoFT2() << " no FT2 entry\n";
      if(ft1.nCutZeroExposure())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutZeroExposure() << " zero exposure\n";
      if(ft1.nCutROI())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutROI() << " outside ROI\n";
      if(ft1.nCutEnergy())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutEnergy() << " outside energy range\n";
    }

  // --------------------------------------------------------------------------
  // COMPUTE PERIODOGRAM
  // --------------------------------------------------------------------------

  if(verbose>=2)
    {
      std::cout << "Computing periodogram...\n";
    }

  Scanner scanner(scan0, scanN, nscan, ft2exp, verbose, log_scan,
		  phase_scan, phase_scan_freq);
  scanner.scan(ft1.events);

  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<IRFs*>::iterator iirf=irfs.begin();iirf!=irfs.end();iirf++)
    delete *iirf;
  for(std::vector<PLIntegratedEA*>::iterator ie = eaint.begin();
      ie!=eaint.end(); ie++)delete *ie;
}
