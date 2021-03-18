//-*-mode:c++; mode:font-lock;-*-

// flareprobability - Stephen Fegan 
//                  - sfegan@llr.in2p3.fr
//                  - 2010-08-03
//
// Calculate probability of photon time separation assuming no flaring
//
// $Id: flare_probability.cpp 4821 2012-12-02 17:36:29Z sfegan $

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

#include "FT1ROIFilter.hpp"
#include "FT2Exp.hpp"
#include "IRF.hpp"
#include "PLIntegratedEA.hpp"
#include "GTI.hpp"
#include "FITS.hpp"
#include "Util.hpp"
#include "Accumulator.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

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

  unsigned nstack=2;
  options.findWithValue("ndifferences", nstack,
			"Number of photon time differences to consider.");
  if(nstack<1)nstack=1;
 
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

  std::string exp_blob_fn = std::string();
  options.findWithValue("exp_file", exp_blob_fn,
			"Load exposure calculation to a file or load it "
			"back from the file if it exists.","exposure");  

  std::string exp_lc_fn = std::string();
  options.findWithValue("exp_lc", exp_lc_fn,
			"Specify the name of a file in which to write the "
			"lightcurve of exposure values. This is mostly "
			"meant as a diagnostics that can be compared with a "
			"curve produced by gtexposure.","exposure");

  double exp_lc_period = 1.0;
  options.findWithValue("exp_lc_period", exp_lc_period,
			"Specify the period of the exposure lightcurve in "
			"days.", "exposure");

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

      if(emax!=0 && energy_cuts.second>emax)
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

  irfLoader::Loader::go();
  const std::vector<std::string> & 
    irf_names(irfLoader::Loader::respIds().find(irf)->second);
  irfInterface::IrfsFactory & factory(*irfInterface::IrfsFactory::instance());
  std::vector<irfInterface::Irfs *> irfs(irf_names.size());
  for (size_t iirf = 0; iirf < irf_names.size(); iirf++)
    irfs[iirf] = factory.create(irf_names[iirf]);

  std::vector<PLIntegratedEA*> eaint(irfs.size());
  for (size_t iirf = 0; iirf < irfs.size(); iirf++)
    eaint[iirf] = new PLIntegratedEA(irfs[iirf], emin, emax, gamma);

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
  
  FT2Exp ft2exp(eaint, gti, ra, dec, exp_lc_period);

  bool loaded_from_blob = false;
  if(!exp_blob_fn.empty())
    {
      BLOBUnserializer blob(exp_blob_fn);
      if(verbose>=2)
	std::cout << "Loading exposure BLOB: " << exp_blob_fn;
      if(blob.good())
	loaded_from_blob = ft2exp.partiallyUnserializeFromBlob(blob);
      if(loaded_from_blob && verbose>=2)std::cout << " .. succeeded\n";
      else if(verbose>=2)std::cout << " .. failed\n";
    }

  if(!loaded_from_blob)
    {
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
	    << "FT2: Total time exposure: " 
	    << lrint(ft2exp.tInt().back().livetime_int_stop)
	    << " seconds\n";
	  for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
	    std::cout
	      << "FT2: Total exposure (EC=" << iea << "): " 
	      << ft2exp.eInt(iea).back().stop << " cm^2 seconds\n";
	}

      if(!exp_blob_fn.empty())
	{
	  BLOBSerializer blob(exp_blob_fn);
	  if(verbose>=2)
	    std::cout << "Saving exposure BLOB: " 
		      << exp_blob_fn << '\n';
	  ft2exp.partiallySerializeToBlob(blob);
	}
    }

  if(!exp_lc_fn.empty())
    {
      const std::vector<std::pair<double,double> >& lc(ft2exp.getExpLC());
      if(verbose>=2)
	std::cout
	  << "FT2: writing exposure lightcurve to: " << exp_lc_fn << '\n';
      std::ofstream str(exp_lc_fn.c_str());
      for(unsigned ilc=0;ilc<lc.size();ilc++)
	str << std::setprecision(10) 
	    << lc[ilc].first << ' ' << lc[ilc].second << '\n';
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
    }

  // --------------------------------------------------------------------------
  // DO CALCULATION
  // --------------------------------------------------------------------------

  double total_exposure = 0;
  for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
    total_exposure += ft2exp.eInt(iea).back().stop;
  double mean_rate = double(ft1.events.size())/total_exposure;
  if(verbose>=1)
    std::cout << "Mean rate: " << mean_rate << " ph/cm^2/s\n";

  std::vector<double> event_stack;
  event_stack.reserve(nstack);
  for(unsigned ievent=0;ievent<ft1.events.size();ievent++)
    {
      unsigned ift2 = ft1.events[ievent].ift2exp;
      double x = 
	(ft1.events[ievent].time-ft2exp.t0()-ft2exp.tInt()[ift2].start)
	/(ft2exp.tInt()[ift2].stop-ft2exp.tInt()[ift2].start);
      double t = 0;
      for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
	t += ft2exp.eInt(iea)[ift2].start + ft2exp.eInt(iea)[ift2].dexp*x;
      std::cout 
	<< std::fixed
	<< lrint(ft1.events[ievent].time) << ' ' 
	<< std::setprecision(3) << ft1.events[ievent].time/86400+51910 << ' '
	<< std::setw(10) << ft1.events[ievent].energy << ' '
	<< std::scientific << std::setprecision(7) << t;
      for(unsigned istack=0;istack<event_stack.size();istack++)
	{
	  double dt = t-event_stack[istack];
	  std::cout << ' ' << std::scientific << std::setprecision(3) << dt;
	  double p = gsl_sf_gamma_inc_P(double(istack+1), mean_rate*dt);
	  std::cout << ' ' << p;
	}
      if(event_stack.size() != nstack)
	event_stack.resize(event_stack.size()+1);
      for(unsigned istack=event_stack.size()-1;istack!=0;istack--)
	event_stack[istack] = event_stack[istack-1];
      event_stack[0] = t;
      std::cout << '\n';
    }

  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<irfInterface::Irfs *>::iterator iirf=irfs.begin();
      iirf!=irfs.end();iirf++)delete *iirf;
  for(std::vector<PLIntegratedEA*>::iterator ie = eaint.begin();
      ie!=eaint.end(); ie++)delete *ie;
}
