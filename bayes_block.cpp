//-*-mode:c++; mode:font-lock;-*-

// bayes_block - Stephen Fegan 
//             - sfegan@llr.in2p3.fr
//             - 2010-08-25
//
// Partition FT1 photons into constant rate blocks
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
#include "PLIntegratedEA.hpp"
#include "GTI.hpp"
#include "FITS.hpp"
#include "Util.hpp"
#include "Accumulator.hpp"
#include "BlockPartition.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

typedef FT2ROI::Interval Interval;

class PoissonLogLikeBlockCost: public BlockCostFunction<Interval>
{
public:
  PoissonLogLikeBlockCost(): BlockCostFunction<Interval>() { }  
  virtual ~PoissonLogLikeBlockCost();
  virtual double cost(std::vector<Interval>::const_iterator begin,
		      std::vector<Interval>::const_iterator end);
};

PoissonLogLikeBlockCost::~PoissonLogLikeBlockCost()
{
  // nothing to see here
}

double PoissonLogLikeBlockCost::
cost(std::vector<Interval>::const_iterator begin,
     std::vector<Interval>::const_iterator end)
{
  const double N = double(end-begin);
  const double M = (end-1)->stop - begin->start;
  return N*(log(M)-log(N));
};


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

void events2intervals(std::vector<Interval>& intervals,
		      const std::vector<FT1AEff> events,
		      const FT2Exp& ft2exp, double mean_rate)
{
  intervals.resize(events.size());
  intervals[0].start=0;
  for(unsigned ievent=0;ievent<events.size();ievent++)
    {
      const double t = events[ievent].time-ft2exp.t0();
      const unsigned iint = events[ievent].ift2exp;
      const FT2Exp::Interval& T(ft2exp.tInt()[iint]);
      const double x = (t-T.start)/(T.stop-T.start);

      double e = 0;
      for(unsigned ie=0;ie<ft2exp.nEInt();ie++)
	if(!ft2exp.eInt(ie).empty())
	  {
	    const FT2Exp::Interval& E(ft2exp.eInt(ie)[iint]);
	    e += E.start + x*(E.stop-E.start);
	  }
      e *= mean_rate;

      double elast = 0;
      if(ievent)elast = intervals[ievent-1].stop;

#if 0
      if(ievent > 17730 && ievent < 17750)
	std::cout
	  << ievent << ' ' << lrint(t) << ' ' 
	  << r2d(events[ievent].theta) << ' ' << events[ievent].energy << ' ' 
	  << events[ievent].event_id << ' ' << events[ievent].run_id << ' ' 
	  << iint << ' '
	  << ft2exp.eInt(0)[iint].stop-ft2exp.eInt(0)[iint].start << ' '
	  << x << ' ' << e << ' ' << e-elast << '\n';
#endif 

      intervals[ievent].start = elast;
      intervals[ievent].stop  = e;
    }
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

  double penalty = 0.0;
  options.findWithValue("penalty", penalty,
			"Block creation penalty.");

  unsigned nsim = 0;
  options.findWithValue("nsim", nsim,
			"Number of simulations to perform.");

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
			      verbose>=2, caldb);
    }

  double eaint_ctheta0 = cos(tmax);
  double eaint_dctheta = 0.0025;

  std::vector<PLIntegratedEA*> eaint(irfs.size());
  for(unsigned iirf=0;iirf<irfs.size();iirf++)
    {
      if(irfs[iirf]->ea())eaint[iirf] = 
	  new PLIntegratedEA(irfs[iirf], eaint_dctheta, eaint_ctheta0,
			     logemin, logemax, gamma, roi);
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
  
  FT2Exp ft2exp(eaint, gti, ra, dec);
  ft2exp.setZMaxCut(zmax);
  ft2exp.setTMaxCut(tmax);
  ft2exp.setTMinCut(tmin);
  ft2exp.setRequireNominalLATMode();

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

  ft2exp.mergeExposures();

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
  std::sort(ft1.events.begin(),ft1.events.end());
  
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
      if(ft1.nCutDuplicate())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutDuplicate() << " duplicated events\n";
    }

  // --------------------------------------------------------------------------
  // CALCULATE EVENT RATES FOR EACH TYPE
  // --------------------------------------------------------------------------

  RandomNumbers* rng = 0;

  double total_exposure = 0;
  for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
    total_exposure += ft2exp.eInt(iea).back().stop;

#if 1
  std::vector<double> rates(ft2exp.nEInt());
  for(unsigned ie=0;ie<ft2exp.nEInt();ie++)
    if(!ft2exp.eInt(ie).empty())
      rates[ie] = ft1.nEvents(ie)/ft2exp.eInt(ie).back().stop;
#else
  std::vector<double> rates(ft2exp.nEInt());
  for(unsigned ie=0;ie<ft2exp.nEInt();ie++)
    if(!ft2exp.eInt(ie).empty())
      rates[ie] = ft1.events.size()/total_exposure;
#endif

  // --------------------------------------------------------------------------
  // SIMULATIONS
  // --------------------------------------------------------------------------

  std::vector<unsigned> sim_npartition;
  std::vector<unsigned> sim_int_npartition;
  VSSimpleStat2<double> sim_stat_nevent;
  VSSimpleStat2<double> sim_stat_cost;
  VSSimpleStat2<double> sim_stat_npartition;
  if(nsim>0)
    {
      if(verbose>=1)
	std::cout << "Doing event simulations: " << nsim << '\n';
      if(!rng)rng = new RandomNumbers();
      SimpleEventTimesSimulation sim(*rng, ft2exp, rates);
      for(unsigned isim=0;isim<nsim;isim++)
	{
	  std::vector<FT1AEff> events;
	  sim.generateEventTimes(events);
	  sim_stat_nevent.accumulate(double(events.size()));
	  double mean_rate = double(events.size())/total_exposure;
	  std::vector<Interval> intervals;
	  events2intervals(intervals, events, ft2exp, mean_rate);
	  PoissonLogLikeBlockCost logl_cost;
	  BlockPenalty<Interval> penalty_cost(logl_cost,penalty);
	  BlockPartionCalculator<Interval> block_calc(intervals, penalty_cost);
	  std::vector<unsigned> partition;
	  double cost;
	  cost = block_calc.compute(partition);
	  const unsigned npartition = partition.size()+1;
	  if(npartition >= sim_npartition.size())
	    sim_npartition.resize(npartition+1);
	  sim_npartition[npartition]++;
	  sim_stat_cost.accumulate(cost);
	  sim_stat_npartition.accumulate(double(npartition));
	}
      
      sim_int_npartition.resize(sim_npartition.size());
      for(unsigned inp=sim_npartition.size(),s=0;inp;inp--)
	{
	  s += sim_npartition[inp-1];
	  sim_int_npartition[inp-1] = s;
	}

      for(unsigned inp=1;inp<sim_npartition.size();inp++)
	std::cout << inp << ' ' << sim_npartition[inp] << ' '
		  << sim_int_npartition[inp] << '\n';
      std::cout
	<< "Mean nevent:     " << sim_stat_nevent.mean() 
	<< " +/- " << sim_stat_nevent.dev() << '\n'
	<< "Mean cost:       " << sim_stat_cost.mean() 
	<< " +/- " << sim_stat_cost.dev() << '\n'
	<< "Mean npartition: " << sim_stat_npartition.mean() 
	<< " +/- " << sim_stat_npartition.dev() << '\n';
  }

  // --------------------------------------------------------------------------
  // CALCULATE OPTIMAL BLOCK PARTITION
  // --------------------------------------------------------------------------

  double mean_rate = double(ft1.events.size())/total_exposure;
  if(verbose>=1)
    std::cout << "Mean rate: " << mean_rate << " ph/cm^2/s\n";

  std::vector<Interval> intervals;
  events2intervals(intervals, ft1.events, ft2exp, mean_rate);

  PoissonLogLikeBlockCost logl_cost;
  BlockPenalty<Interval> penalty_cost(logl_cost,penalty);
  BlockPartionCalculator<Interval> block_calc(intervals, penalty_cost);

  std::vector<unsigned> partition;
  double cost;
  cost = block_calc.compute(partition);
  
  double bstart = 0;
  double tstart = ft2exp.t0()/86400+51910;
  unsigned nstart = 0;
  for(unsigned ip=0;ip<=partition.size();ip++)
    {
      unsigned nend;
      double tend;
      if(ip==partition.size())
	{
	  nend = ft1.events.size();
	  tend = ft1.events.back().time/86400+51910;
	}
      else
	{
	  nend = partition[ip];
	  tend = ft1.events[partition[ip]-1].time/86400+51910;
	}

      const double e = (intervals[nend-1].stop-bstart)/mean_rate;
      std::cout 
	<< ip << ' ' 
	<< nend << ' ' << nend-nstart << ' '
	<< tstart << ' ' << tend << ' ' << tend-tstart << ' '
	<< e << ' ' << double(nend-nstart)/e << '\n';
      bstart = intervals[nend-1].stop;
      tstart = tend;
      nstart = nend;
    }

  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<IRFs*>::iterator iirf=irfs.begin();iirf!=irfs.end();iirf++)
    delete *iirf;
  for(std::vector<PLIntegratedEA*>::iterator ie = eaint.begin();
      ie!=eaint.end(); ie++)delete *ie;
  delete rng;
}
