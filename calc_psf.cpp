//-*-mode:c++; mode:font-lock;-*-

//
// calc_psf.cpp - Stephen Fegan 
//              - sfegan@llr.in2p3.fr
//              - 2009-09-09
//
// Fit PSF curves to stacked data
//
// $Id: calc_psf.cpp 2064 2010-08-24 12:06:58Z sfegan $
//

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"
#include "VERITAS/VSDataConverter.hpp"

#include "FT1.hpp"
#include "FT2.hpp"
#include "FT2ROI.hpp"
#include "Catalog.hpp"
#include "Util.hpp"
#include "MyMinuit.hpp"
#include "PSF.hpp"
#include "Analysis.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

class SourceFT2GTIFilter: public FoundSourceVisitor
{
public:
  SourceFT2GTIFilter(FoundSourceVisitor* visitor, 
		   const std::vector<FT2ROI*>& ft2rois);
  virtual ~SourceFT2GTIFilter();
  virtual void visitFSEvent(unsigned irow, FT1& event,
			    unsigned isrc, double ra, double dec, double dist);
  unsigned nFSEventsDispatched() const { return m_nfsevents; }
private:
  FoundSourceVisitor* m_visitor;
  const std::vector<FT2ROI*>& m_ft2rois;
  unsigned m_nfsevents;
};

SourceFT2GTIFilter::
SourceFT2GTIFilter(FoundSourceVisitor* visitor, 
		 const std::vector<FT2ROI*>& ft2rois):
  FoundSourceVisitor(), 
  m_visitor(visitor), m_ft2rois(ft2rois), m_nfsevents() { }

SourceFT2GTIFilter::~SourceFT2GTIFilter()
{
  // nothing to see here
}

void SourceFT2GTIFilter::
visitFSEvent(unsigned irow, FT1& event,
	     unsigned isrc, double ra, double dec, double dist)
{
  if(isrc<m_ft2rois.size() && m_ft2rois[isrc]->isInFT2GTI(event.time))
    {
      m_visitor->visitFSEvent(irow, event, isrc, ra, dec, dist);
      m_nfsevents++;
    }
#if 0
  else
    std::cout << r2d(event.zenith) << ' ' << r2d(event.theta) << ' '
  	      << r2d(dist) << '\n';
#endif
}

class DistVecMaker: public FoundSourceVisitor
{
public:
  DistVecMaker(): FoundSourceVisitor(), m_d() { }
  virtual ~DistVecMaker();
  virtual void visitFileVector(const std::string& filename,
			       const std::string& tablename,
			       unsigned nrow, FITSHeader& header);
  virtual void visitFSEvent(unsigned irow, FT1& event,
			    unsigned isrc, double ra, double dec, double dist);
  std::vector<double> m_d;
};

DistVecMaker::~DistVecMaker()
{
  // nothing to see here
}

void DistVecMaker::visitFileVector(const std::string& filename,
				   const std::string& tablename,
				   unsigned nrow, FITSHeader& header)
{
  m_d.reserve(m_d.size() + nrow);
}

void DistVecMaker::visitFSEvent(unsigned irow, FT1& event, unsigned isrc, 
				double ra, double dec, double dist)
{
  m_d.push_back(r2d(dist));
}

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] ft1_file catalog_file"
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

  bool interactive = false;
  if(options.find("interactive", "Run MINUIT interactively.")
     !=VSOptions::FS_NOT_FOUND)
    interactive = true;

  std::string interactive_commands;
  if(options.findWithValue("icmd", interactive_commands,
			   "Run MINUIT interactively, with given "
			   "commands.")
     !=VSOptions::FS_NOT_FOUND)
    interactive = true;

  std::string ft2_fn="";
  options.findWithValue("ft2", ft2_fn,
			"Load FT2 file for ROI based zenith cuts.");
  
  Catalog::Filter catalog_cuts;
  options.addCatagory("cat","Source catalog selection:");

  if(options.find("cat_default",
		  "Apply default cuts to catalog sources","cat")
     != VSOptions::FS_NOT_FOUND)
    {
      setCut(catalog_cuts.glat_min, 15.0);
      setCut(catalog_cuts.ts_min, 9.0);
      setCut(catalog_cuts.nassoc_min, 1U);
      setCut(catalog_cuts.max_inter_assoc_dist_max, 0.1);
      setCut(catalog_cuts.min_neighbour_dist_min, 3.0);
    }

  if(options.findWithValue("cat_glat", catalog_cuts.glat_min.second,
			   "Set minimum galactic latitude for catalog "
			   "sources.",
			   "cat")
     != VSOptions::FS_NOT_FOUND)
    catalog_cuts.glat_min.first = true;

  if(options.findWithValue("cat_ts", catalog_cuts.ts_min.second,
			   "Set minimum TS value for catalog sources.",
			   "cat")
     != VSOptions::FS_NOT_FOUND)catalog_cuts.ts_min.first = true;

  options.findWithValue("cat_band_ts", catalog_cuts.band_ts_min,
			"Set minimum band TS values for catalog sources.",
			"cat");

  if(options.findWithValue("cat_nassoc", catalog_cuts.nassoc_min.second,
			   "Set minimum number of associations for "
			   "catalog sources",
			   "cat")
     != VSOptions::FS_NOT_FOUND)catalog_cuts.nassoc_min.first = true;
  
  if(options.findWithValue("cat_assoc_dist", 
			   catalog_cuts.max_inter_assoc_dist_max.second,
			   "Set maximum distance between most widely "
			   "separated assoications for catalog sources. This "
			   "cut can be used to eliminate sources which have "
			   "multiple associations.",
			   "cat")
     != VSOptions::FS_NOT_FOUND)
    catalog_cuts.max_inter_assoc_dist_max.first = true;

  if(options.findWithValue("cat_nbr_dist", 
			   catalog_cuts.min_neighbour_dist_min.second,
			   "Set minimum distance between two catalog sources "
			   "that are above a certain TS threshold (see next "
			   "option).",
			   "cat")
     != VSOptions::FS_NOT_FOUND)
    catalog_cuts.min_neighbour_dist_min.first = true;

  double cat_nbr_ts = 0;
  options.findWithValue("cat_nbr_ts", cat_nbr_ts,
			"Set minimum TS value for neighbours to be cut with "
			"distance cut.", "cat");

  catalog_cuts.glat_min.second = d2r(catalog_cuts.glat_min.second);
  catalog_cuts.max_inter_assoc_dist_max.second =
    d2r(catalog_cuts.max_inter_assoc_dist_max.second);
  catalog_cuts.min_neighbour_dist_min.second =
    d2r(catalog_cuts.min_neighbour_dist_min.second);

  bool catalog_strike_overlapping = true;
  if(options.find("cat_no_strike_overlapping",
		  "Do not delete catalog sources with overlapping ROIs.",
		  "cat")
     != VSOptions::FS_NOT_FOUND)catalog_strike_overlapping = false;
  
  std::set<unsigned> catalog_strike_num;
  std::set<std::string> catalog_strike_name;
  options.findWithValue("cat_strike_num", 
			catalog_strike_num,
			"Strike entries from catalog by number in list (use "
			"dump_catalog option initially to see which source "
			"this corresponds to). Numbers should be given as a "
			"comma separated list.",
			"cat");
  options.findWithValue("cat_strike_name", 
			catalog_strike_name,
			"Strike entries from catalog by name (use "
			"dump_catalog option initially to see which sources "
			"are in the catalog). Names should be given as a "
			"comma separated list.",
			"cat");  
 
  bool dump_catalog = false;
  bool dump_catalog_assoc = false;
  std::string dump_catalog_fn;
  if(options.findWithValue("dump_catalog", dump_catalog_fn,
			   "Dump filtered catalog to stdout or to a file.",
			   "cat")
     != VSOptions::FS_NOT_FOUND)dump_catalog = true;
  if(options.find("dump_assoc", "Dump catalog associations also.",
		  "cat")
     != VSOptions::FS_NOT_FOUND)dump_catalog_assoc = true;

  bool use_catalog_position = false;
  if(options.find("use_catalog_position",
		  "Use the fitted catalog position rather than the "
		  "position of the most probable association.",
		  "cat")
     != VSOptions::FS_NOT_FOUND)use_catalog_position = true;

  if(!use_catalog_position && (catalog_cuts.nassoc_min.first==false ||
			       catalog_cuts.nassoc_min.second==0))
    {
      std::cerr 
	<< "WARNING: The cut on number of associations was automatically set to 1. Use\n"
	<< "WARNING: the \"-use_catalog_position\" option or specify a different value of\n"
	<< "WARNING: this cut using \"-cat_nassoc\".\n";
      catalog_cuts.nassoc_min.first = true;
      catalog_cuts.nassoc_min.second = 1;
    }

  FT1::Filter ft1_cuts;
  options.addCatagory("ft1","Data selection:");
  
  if(options.findWithValue("emin", ft1_cuts.e_min.second,
			   "Set minimum energy cut for events.",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.e_min.first = true;

  if(options.findWithValue("emax", ft1_cuts.e_max.second,
			   "Set maximum energy cut for events.",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.e_max.first = true;

  if(options.findWithValue("logemin", ft1_cuts.e_min.second,
			   "Set minimum energy cut for events (set in log10).",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.e_min.first = true,
      ft1_cuts.e_min.second = pow(10.0,ft1_cuts.e_min.second);
  
  if(options.findWithValue("logemax", ft1_cuts.e_max.second,
			   "Set maximum energy cut for events (set in log10).",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.e_max.first = true,
      ft1_cuts.e_max.second = pow(10.0,ft1_cuts.e_max.second);

  if(options.findWithValue("zmax", ft1_cuts.z_max.second,
			   "Set maximum zenith angle cut for events.",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.z_max.first = true;
  ft1_cuts.z_max.second = d2r(ft1_cuts.z_max.second);

  if(options.findWithValue("tmin", ft1_cuts.t_min.second,
			   "Set minimum theta cut for events.",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.t_min.first = true,
      ft1_cuts.t_min.second = d2r(ft1_cuts.t_min.second);

  if(options.findWithValue("tmax", ft1_cuts.t_max.second,
			   "Set maximum theta cut for events.",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.t_max.first = true,
      ft1_cuts.t_max.second = d2r(ft1_cuts.t_max.second);

  if(options.findWithValue("costmax", ft1_cuts.t_min.second,
			   "Set minimum theta cut for events (set in cosine).",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.t_min.first = true,
      ft1_cuts.t_min.second = std::acos(ft1_cuts.t_min.second);
  
  if(options.findWithValue("costmin", ft1_cuts.t_max.second,
			   "Set maximum theta cut for events (set in cosine).",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.t_max.first = true,
      ft1_cuts.t_max.second = std::acos(ft1_cuts.t_max.second);

  if(options.findWithValue("ev_class", ft1_cuts.event_class.second,
			   "Set required event class.",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.event_class.first = true;
  
  if(options.findWithValue("conv_type", ft1_cuts.conversion_type.second,
			   "Set required conversion type (0=front, 1=back).",
			   "ft1")
     != VSOptions::FS_NOT_FOUND)
    ft1_cuts.conversion_type.first = true;

  options.addCatagory("catfit","Catalog source matching:");
  double d_max = 1.5;
  options.findWithValue("dmax", d_max,
			"Set maximum distance that an event can be from "
			"a catalog source and still be consideted as "
			"coming from that source.","catfit");
  d_max = d2r(d_max);

  options.addCatagory("psffit","Functional forms of the PSF:");
  bool dblking_psf_fit = false;
  if(options.find("dblking_psf",
		  "Fit a double King function to the integral PSF "
		  "distribution (rather than a single King).",
		  "psffit")
     != VSOptions::FS_NOT_FOUND)dblking_psf_fit = true;  

  bool gaussian_psf_fit = false;
  if(options.find("gaussian_psf",
		  "Fit a Gaussian function to the integral PSF "
		  "distribution (rather than a King).",
		  "psffit")
     != VSOptions::FS_NOT_FOUND)gaussian_psf_fit = true;  

  bool dblgaussian_psf_fit = false;
  if(options.find("dblgaussian_psf",
		  "Fit a double Gaussian function to the integral "
		  "PSF distribution (rather than a King).",
		  "psffit")
     != VSOptions::FS_NOT_FOUND)dblgaussian_psf_fit = true;  

  unsigned multigaussian_n = 0;
  bool multigaussian_psf_fit = false;
  if(options.findWithValue("multigaussian_psf", multigaussian_n,
			   "Fit a multiple Gaussian function to the integral "
			   "PSF distribution (rather than a King). "
			   "Specify the number of Gaussian functions to use.",
			   "psffit")
     != VSOptions::FS_NOT_FOUND)multigaussian_psf_fit = true;  

  bool dump_psf_fit = false;
  std::string dump_psf_fit_fn = "";
  if(options.findWithValue("dump_psf", dump_psf_fit_fn,
			   "Dump integral PSF profile from data and fit "
			   "to stdout or to a file.",
			   "psffit")
     != VSOptions::FS_NOT_FOUND)dump_psf_fit = true;  

  std::string eodat_fn="";
  options.findWithValue("eodat", eodat_fn,
			"Load event orientation data file and refit event "
			"direction. Use the \"calc_event_orient\" program "
			"to generate this data file.",
			"refit");

  triple<double,double,double> sc_rot_triple;
  options.findWithValue("scrot", sc_rot_triple,
			"Specify event orientation rotation. Should be given "
			"as three numbers \"rx/ry/rz\", with each number in "
			"arcminutes.",
			"refit");
  Vec3D sc_rot(0,0,0);
  sc_rot &= Vec3D(0,0,d2r(sc_rot_triple.third/60.0));
  sc_rot &= Vec3D(0,d2r(sc_rot_triple.second/60.0),0);
  sc_rot &= Vec3D(d2r(sc_rot_triple.first/60.0),0,0);

  std::vector<double> sc_fourier_theta;
  options.findWithValue("theta_coeff", sc_fourier_theta,
			"Specify coeffients of (odd) Fourier series expansion "
			"of theta correction. Coefficients should be given "
			"in arcminutes.",
			"refit");
  for(unsigned ift=0;ift<sc_fourier_theta.size();ift++)
    sc_fourier_theta[ift] = d2r(sc_fourier_theta[ift]);

  bool do_radec_rot_triple = false;
  triple<double,double,double> radec_rot_triple;
  Vec3D radec_rot(0,0,0);
  if(options.findWithValue("radecrot", radec_rot_triple,
			   "Specify RADEC rotation. Should be given as"
			   "three numbers \"rx/ry/rz\", with each number in "
			   "arcminutes.",
			   "refit")
     != VSOptions::FS_NOT_FOUND)
    {
      do_radec_rot_triple = true;
      radec_rot &= Vec3D(0,0,d2r(radec_rot_triple.third/60.0));
      radec_rot &= Vec3D(0,d2r(radec_rot_triple.second/60.0),0);
      radec_rot &= Vec3D(d2r(radec_rot_triple.first/60.0),0,0);
    }
  
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

  int arg_req = 2;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string ft1_fn(*argv);
  argc--, argv++;

  std::string catalog_fn(*argv);
  argc--, argv++;

  // --------------------------------------------------------------------------
  // LOAD AND FILTER CATALOG
  // --------------------------------------------------------------------------

  if(verbose>0)std::cout << command_line << '\n';

  std::vector<Catalog> catalog_vec;
  unsigned ncat = Catalog::loadFromFITS(catalog_vec, catalog_fn,
					Catalog::tableName(), cat_nbr_ts);
  if(verbose>0)
    std::cout << "Loaded catalog (" << catalog_fn << ") has " 
	      << ncat << " entries\n";

  if(catalog_cuts.haveCuts())
    {
      if(verbose>0)
	{
	  std::cout << "Applying cuts to catalog: \n";
	  catalog_cuts.print(std::cout);
	}
      
      ncat = Catalog::filterCatalog(catalog_vec, catalog_cuts);

      if(verbose>0)
	std::cout << "Filtered catalog has " << ncat << " entries\n";
    }
  else if(verbose>0)
    std::cout << "Loaded catalog has " << ncat << " entries\n";

  unsigned oldncat = ncat;

  std::set<unsigned> ols_catalog_strike_num;

  for(unsigned icat=0;icat<ncat;icat++)
    for(unsigned jcat=icat+1;jcat<ncat;jcat++)
      {
	double d = 
	  sphere_dist(catalog_vec[icat].ra, catalog_vec[icat].dec,
		      catalog_vec[jcat].ra, catalog_vec[jcat].dec);
	if(d<2.0*d_max)
	  {
	    if(!catalog_strike_overlapping)
	      std::cout 
		<< "*Warning* Overlapping ROIs: " 
		<< catalog_vec[icat].nickname << " [" << icat+1 << "] and " 
		<< catalog_vec[jcat].nickname << " [" << jcat+1 << "] with d=" 
		<< r2d(d) << " deg\n";
	    else
	      {
		ols_catalog_strike_num.insert(icat);
		ols_catalog_strike_num.insert(jcat);
	      }
	  }
      }

  for(std::set<unsigned>::const_reverse_iterator istrike = 
	ols_catalog_strike_num.rbegin();
      istrike != ols_catalog_strike_num.rend(); istrike++)
    {
      unsigned strike = *istrike;
      if(strike>0 && strike<=ncat)
	{
	  strike--;
	  catalog_vec.erase(catalog_vec.begin()+strike);
	  ncat--;
	}
    }

  if(verbose>0 && oldncat!=ncat)
    std::cout << "Non-overlapping catalog has " << ncat << " entries\n";

  oldncat=ncat;

  if(dump_catalog)
    {
      if(dump_catalog_fn.empty() 
	 || dump_catalog_fn=="-" || dump_catalog_fn=="STDOUT")
	Catalog::printCatalog(catalog_vec, std::cout, dump_catalog_assoc);
      else
	{
	  std::ofstream str(dump_catalog_fn.c_str());
	  Catalog::printCatalog(catalog_vec, str, dump_catalog_assoc);
	}
    }

  for(std::set<unsigned>::const_reverse_iterator istrike = 
	catalog_strike_num.rbegin();
      istrike != catalog_strike_num.rend(); istrike++)
    {
      unsigned strike = *istrike;
      if(strike>0 && strike<=ncat)
	{
	  strike--;
	  if(verbose>0)
	    std::cout << "Deleting catalog source: "
		      << catalog_vec[strike].nickname << '\n';
	  catalog_vec.erase(catalog_vec.begin()+strike);
	  ncat--;
	}
    }

  for(std::set<std::string>::const_iterator istrike = 
	catalog_strike_name.begin();
      istrike != catalog_strike_name.end(); istrike++)
    {
      const std::string& name = *istrike;
      std::vector<Catalog>::iterator icat = catalog_vec.begin();
      while(icat != catalog_vec.end() && icat->nickname != name)icat++;
      if(icat != catalog_vec.end())
	{
	  if(verbose>0)
	    std::cout << "Deleting catalog source: " << name << '\n';
	  catalog_vec.erase(icat);
	  ncat--;
	}
    }

  if(oldncat!=ncat)
    {
      if(verbose >= 0)
	std::cout << "Final catalog has " << ncat << " entries\n";

      if(dump_catalog)
	{
	  if(dump_catalog_fn=="-" || dump_catalog_fn=="STDOUT")
	    Catalog::printCatalog(catalog_vec, std::cout, dump_catalog_assoc);
	  else
	    {
	      std::ofstream str(dump_catalog_fn.c_str());
	      Catalog::printCatalog(catalog_vec, str, dump_catalog_assoc);
	    }
	}
    }

  std::vector<double> cat_ra(ncat);
  std::vector<double> cat_dec(ncat);

  if(use_catalog_position)
    for(unsigned icat=0;icat<ncat;icat++)
      cat_ra[icat] = catalog_vec[icat].ra,
	cat_dec[icat] = catalog_vec[icat].dec;
  else
    for(unsigned icat=0;icat<ncat;icat++)
      cat_ra[icat] = catalog_vec[icat].assoc[0].ra,
	cat_dec[icat] = catalog_vec[icat].assoc[0].dec;  

  // --------------------------------------------------------------------------
  // READ IN THE GTI AND LOAD THE SPACECRAFT (FT2) FILE
  // --------------------------------------------------------------------------

  GTIRange gti;
  std::vector<FT2ROI*> ft2rois;
  if((ft1_cuts.z_max.first || ft1_cuts.t_max.first || ft1_cuts.t_min.first)
     && !ft2_fn.empty())
    {
      if(verbose>0)
	std::cout << "Load: GTI - " << ft1_fn << '\n';

      gti.loadGTIsFromFITS(ft1_fn);
      
      if(verbose>0)
	{
	  std::cout << "GTI: Total time " << lrint(gti.totalTime()) 
		    << " seconds in " << gti.nGTIs() << " periods\n";

	  std::cout << "Load: FT2 - " << ft2_fn << '\n';

	  std::cout
	    << "Applying FT2 based cuts, requiring centers of ROIs to have:\n";
	  if(ft1_cuts.z_max.first)
	    std::cout << "- Zenith               <= " 
		      << r2d(ft1_cuts.z_max.second-d_max) << '\n';
	  if(ft1_cuts.t_max.first)
	    std::cout << "- Theta                <= "
		      << r2d(ft1_cuts.t_max.second) << '\n';
	  if(ft1_cuts.t_min.first)
	    std::cout << "- Theta                >= " 
		      << r2d(ft1_cuts.t_min.second) << '\n';
	}

      ft2rois.resize(ncat);
      MultiFITSVectorVisitor<FT2> mv;
      for(unsigned icat=0;icat<ncat;icat++)
	{
	  ft2rois[icat] = new FT2ROI(gti, cat_ra[icat], cat_dec[icat]);
	  ft2rois[icat]->setRequireNominalLATMode();
	  if(ft1_cuts.z_max.first)
	    ft2rois[icat]->setZMaxCut(ft1_cuts.z_max.second-d_max);
	  if(ft1_cuts.t_max.first)
	    ft2rois[icat]->setTMaxCut(ft1_cuts.t_max.second);
	  if(ft1_cuts.t_min.first)
	    ft2rois[icat]->setTMinCut(ft1_cuts.t_min.second);
	  mv.addVisitor(ft2rois[icat]);
	}

      // Zenith and theta cuts are handled by FT2 not by FT1 - but can
      // use FT1 cuts to speed things up hand by cutting events early
      if(ft1_cuts.t_max.first)
	ft1_cuts.t_max.second += d_max;
      if(ft1_cuts.t_min.first && ft1_cuts.t_min.second>d_max)
	ft1_cuts.t_min.second -= d_max;
      else
	ft1_cuts.t_min.first = false;

      FITSVectorDispatcher<FT2> dispatcher(&mv);
      dispatcher.dispatchVector(ft2_fn, FT2::tableName());
      if(verbose>0)
	for(unsigned icat=0;icat<ncat;icat++)
	  {
	    std::cout << "- ROI exposure: " << catalog_vec[icat].nickname 
		      << " = " << lrint(ft2rois[icat]->gti().totalTime())
		      << '\n';
	  }
    }

  // --------------------------------------------------------------------------
  // PROCESS THE EVENTS
  // --------------------------------------------------------------------------

  FITSVectorVisitor<FT1>* visitor(0);

  DistVecMaker* dvm(new DistVecMaker);
  FoundSourceVisitor* fsv(dvm);
  SourceFT2GTIFilter* ft2filter(0);
  if(!ft2rois.empty())
    {
      ft2filter = new SourceFT2GTIFilter(dvm, ft2rois);
      fsv = ft2filter;
    }
  SourceFinder* finder(new SourceFinder(fsv, d_max, cat_ra, cat_dec));
  visitor = finder;

  RADECRotation* radec_rotator(0);
  if(do_radec_rot_triple)
    {
      radec_rotator = new RADECRotation(visitor, radec_rot);
      visitor = radec_rotator;
    }

  EventRADECRecalc* radec_recalc(0);
  EODatEntryMatcher* eodat_matcher(0);
  if(!eodat_fn.empty())
    {
      radec_recalc = new EventRADECRecalc(visitor, sc_rot, sc_fourier_theta);
      eodat_matcher = new EODatEntryMatcher(radec_recalc, eodat_fn);
      visitor = eodat_matcher;
    }

  FT1FilterVisitor* ft1_filter(0);
  if(ft1_cuts.haveCuts())
    {
      if(verbose>0)
	{
	  std::cout << "Applying cuts to events: \n";
	  ft1_cuts.print(std::cout);
	}
      ft1_filter = new FT1FilterVisitor(visitor, ft1_cuts);
      visitor = ft1_filter;
    }

  FITSVectorDispatcher<FT1> dispatcher(visitor);
  unsigned nraw = 
    dispatcher.dispatchVector(ft1_fn, FT1::tableName());

  if(verbose>0)
    {
      std::cout << "Processed: " << nraw << " raw events\n";
      if(ft1_filter)
	std::cout << "           " << ft1_filter->nElementsDispatched()
		  << " after selection cuts\n";
      if(eodat_matcher)
	std::cout << "           " << eodat_matcher->nMatchedEventsDispatched()
		  << " after event orientation matching\n";
	
      std::cout << "           " << finder->nFSEventsDispatched() 
		<< " with valid catalog source\n";

      if(ft2filter)
	std::cout << "           " << ft2filter->nFSEventsDispatched() 
		  << " with valid ROI FT2 interval\n";
    }

  std::vector<double> d_vec = dvm->m_d;

  std::sort(d_vec.begin(),d_vec.end());
  unsigned nd = d_vec.size();

  PSFCalc* psf_calc(0);
  if(multigaussian_psf_fit)
    psf_calc = new Like_MultiGaussianPSFWithBkg(multigaussian_n,
						d_vec, r2d(d_max));
  else if(dblgaussian_psf_fit)
    psf_calc = new Like_DblGaussianPSFWithBkg(d_vec, r2d(d_max));
  else if(gaussian_psf_fit)
    psf_calc = new Like_GaussianPSFWithBkg(d_vec, r2d(d_max));
  else if(dblking_psf_fit)
    psf_calc = new Like_DblKingPSFWithBkg(d_vec, r2d(d_max));
  else
    //psf_calc = new Like_BurdettPSFWithBkg(d_vec, r2d(d_max));
    psf_calc = new Like_KingPSFWithBkg(d_vec, r2d(d_max));
  
  MyMinuit minuit(psf_calc, verbose-2);
  minuit.minimize(interactive, 1e-3, 0.5, interactive_commands);
  MyMinuit::MinStat ms = minuit.getMinimizationStatus();

  unsigned np = psf_calc->numParam();
  std::vector<double> p = minuit.pVals();
  std::vector<double> p_err = minuit.pErrs();
  std::vector<std::vector<double> > p_cov = minuit.pCovMtx();
  
  double r68, r68_err;
  psf_calc->modelRP(r68,r68_err,p,p_err,p_cov,0.68);

  PSFCalc* bkg_fn = new Like_IsotropicBkg(d_vec, r2d(d_max));
  std::vector<double> bkg_p(1,1.0);

  std::vector<double> x(1);
  std::vector<double> x0(1,0);
  std::vector<double> P_vec(nd);
  std::vector<double> bkg_P_vec(nd);
  for(unsigned id=0;id<nd;id++)
    {
      x[0] = d_vec[id];
      P_vec[id] = psf_calc->Pint(x,x0,p);
      bkg_P_vec[id] = bkg_fn->Pint(x,x0,bkg_p);
    }

  if(dump_psf_fit)
    {
      std::ostream* str(0);
      std::ofstream* fstr(0);

      if(dump_psf_fit_fn.empty() 
	 || dump_psf_fit_fn=="-" || dump_psf_fit_fn=="STDOUT")	
	str = &std::cout;
      else
	str = fstr = new std::ofstream(dump_psf_fit_fn.c_str());

      for(unsigned id=0;id<nd;id++)
	(*str) << d_vec[id] << ' ' << id << ' ' 
	       << P_vec[id] << ' ' << bkg_P_vec[id] << '\n';
      
      delete fstr;
    }

  if(verbose>0)
    {
      for(unsigned ip=0;ip<np;ip++)
	{
	  std::cout << p[ip] << '\t' << p_err[ip];
	  for(unsigned jp=0;jp<np;jp++)
	    std::cout << '\t' << p_cov[ip][jp];
	  std::cout << '\n';
	}
      
      std::cout << "r68 = " << r68 << " +/- " << r68_err
		<< " LogL = " << ms.fmin << ' ' << ms.istat 
		<< " (bkg only LogL = " << bkg_fn->f(bkg_p) << ")\n";
    }
  else
    {
      std::cout << VSDataConverter::toString(r68) << ' '
		<< VSDataConverter::toString(r68_err) << '\n';
    }
}
