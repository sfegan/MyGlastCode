//-*-mode:c++; mode:font-lock;-*-

// exposure_comparison - Stephen Fegan 
//                     - sfegan@llr.in2p3.fr
//                     - 2011-05-18
//
// Compare ST and my exposure calculation
//
// $Id: flare_probability.cpp 2034 2010-08-04 12:04:55Z sfegan $

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"
#include "VERITAS/VSFileUtility.hpp"
#include "VERITAS/VSDataConverter.hpp"
#include "VERITAS/VSSimpleStat.hpp"

#include "FT1ROIFilter.hpp"
#include "FT2Exp.hpp"
#include "IRF.hpp"
#include "GTI.hpp"
#include "FITS.hpp"
#include "Util.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

typedef FT2ROI::Interval Interval;

class FT2ExpComp: public FT2ROI
{
public:

  FT2ExpComp(const std::vector<IRFs*>& irfs,
	     const GTIRange& gti, double ra, double dec): 
    FT2ROI(gti,ra,dec), m_irfs(irfs)
  {
    setRequireNominalLATMode();
  }
  ~FT2ExpComp();

  virtual void visitAcceptedInterval(unsigned irow, FT2& ft2);

private:
  const std::vector<IRFs*>&           m_irfs;
};

void FT2ExpComp::visitAcceptedInterval(unsigned irow, FT2& ft2)
{
  double dt = ft2.livetime;
  double livetime_frac = dt/(ft2.t_stop-ft2.t_start);
  m_tacc.add(dt);
  m_tint.push_back(Interval(ft2.t_start-m_t0,ft2.t_stop-m_t0,dt,m_tacc.sum()));

  const double ct = 
    std::cos(sphere_dist(m_ra, m_dec, ft2.scz_ra, ft2.scz_dec));

  unsigned day = lrint(ft2.t_start/86400.0-0.5);
  if(day != m_day)
    {
      if(m_day != 0)m_dayexp.push_back(std::make_pair<unsigned,double>
				       (m_day*86400+43200,m_dayacc.sum()));
      m_dayacc.reset();
      m_day = day;
    }

  for(unsigned iea=0;iea<m_ea.size();iea++)
    if(m_ea[iea])
      {
	double aeff = m_ea[iea]->value(ct, 0, livetime_frac)*10000;

	Interval i;
	i.start = m_eacc[iea].sum();
	i.x1 = aeff;
	double exposure = aeff*dt;
	i.x2 = exposure;
	m_eacc[iea].add(exposure);
	i.stop = m_eacc[iea].sum();
	m_eint[iea].push_back(i);
	m_dayacc.add(exposure);
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
  if(options.find("v","Print verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 1;
  if(options.find("vv","Print very verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 2;

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

  double zmax = 105;
  options.findWithValue("zmax", zmax,
			"Set maximum zenith angle cut for events.",
			"cuts");
  zmax = d2r(zmax);

  double tmin = 0;
  options.findWithValue("tmin", tmin,
			"Set minimum theta cut for events.",
			"cuts");
  tmin = d2r(tmin);

  double tmax = 60;
  options.findWithValue("tmax", tmax,
			"Set maximum theta cut for events.",
			"cuts");
  tmax = d2r(tmax);

  // -------------------------------- EXPOSURE --------------------------------

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
  if(verbose>=2)
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

  if(verbose>=2)
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
			      verbose>=2, caldb);
    }

  // --------------------------------------------------------------------------
  // READ IN FT1 GTI
  // --------------------------------------------------------------------------

  GTIRange gti;
  if(verbose>=2)
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
  
  FT2ExpComp ft2exp(irfs, gti, ra, dec);
  ft2exp.setZMaxCut(zmax);
  ft2exp.setTMaxCut(tmax);
  ft2exp.setTMinCut(tmin);
  ft2exp.setRequireNominalLATMode();

  FITSVectorDispatcher<FT2> ft2_dispatcher(&ft2exp);

  if(verbose>=2)
    std::cout << "Load: FT2 - " << ft2_fn << '\n';
  ft2_dispatcher.dispatchVector(ft2_fn, FT2::tableName());

  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<IRFs*>::iterator iirf=irfs.begin();iirf!=irfs.end();iirf++)
    delete *iirf;
}
