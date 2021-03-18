//-*-mode:c++; mode:font-lock;-*-

// avg_psf.cpp - Stephen Fegan 
//             - sfegan@llr.in2p3.fr
//             - 2010-07-28
//
// Average the PSF vs ENERGY over THETA angles using the Aeff and PSF IRFs
//
// $Id: avg_psf.cpp 1995 2010-07-26 14:33:48Z sfegan $

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSFileUtility.hpp"

#include "IRF.hpp"
#include "Util.hpp"
#include "Accumulator.hpp"

using namespace VERITAS;

// ****************************************************************************
//
// MAIN
//
// ****************************************************************************

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname << " irfs"
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
  if(options.find("q","Print no messages during running.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 0;
  if(options.find("v","Print verbose messages during running.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 1;

  double tmin = 0;
  options.findWithValue("tmin", tmin,
			"Set desired minimum theta angle.");

  double tmax = 60;
  options.findWithValue("tmax", tmax,
			"Set desired maximum theta angle.");

  const double costmin = std::cos(d2r(tmin));
  const double costmax = std::cos(d2r(tmax));

  std::string caldb = "$CALDB";
  options.findWithValue("caldb", caldb,
			"Base directory of calibration database.");

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

  int arg_req = 1;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string irf(*argv);
  argc--, argv++;  

  // --------------------------------------------------------------------------
  // READ IN IRFs
  // --------------------------------------------------------------------------

  std::vector<IRFs*> irfs(2);
  for(unsigned iirf=0;iirf<irfs.size();iirf++)
    {
      irfs[iirf] = new IRFs;
      irfs[iirf]->loadAllIRFs(irf, iirf==0?true:false, true, verbose>=1, caldb);
    }

  // --------------------------------------------------------------------------
  // DO CALCULATIONS
  // --------------------------------------------------------------------------

  const double dd = 0.01;
  std::vector<double> d;
  std::vector<double> dA;
  for(double deg=0.5*dd;deg<10;deg+=dd)
    {
      const double rdi = d2r(deg-0.5*dd);
      const double rdo = d2r(deg+0.5*dd);      
      d.push_back(d2r(deg));
      dA.push_back(M_PI*(rdo*rdo - rdi*rdi));
    }

  const double dct = (costmin-costmax)*0.01;
  std::vector<double> ct;
  for(double costheta=costmax+0.5*dct;costheta<costmin;costheta+=dct)
    ct.push_back(costheta);

  for(double log10e=1;log10e<6.001;log10e+=0.05)
    {
      std::cout << log10e << ' ' << std::pow(10.0,log10e);
      for(unsigned iirf=0;iirf<irfs.size();iirf++)
	{
	  std::vector<double> ea(ct.size());
	  //	  Accumulator aea;
	  for(unsigned ict=0;ict<ct.size();ict++)
	    {
	      const double a = irfs[iirf]->ea()->value(log10e, ct[ict]);
	      ea[ict]=a;
	      //	      aea.add(a);
	    }

	  std::vector<double> avgpsf;
	  Accumulator aavgpsf;
	  for(unsigned id=0;id<d.size();id++)
	    {
	      Accumulator ap;
	      for(unsigned ict=0;ict<ct.size();ict++)
		{
		  const double p = 
		    irfs[iirf]->psf()->value(d[id],log10e,ct[ict]);
		  ap.add(p*ea[ict]);
		}
	      aavgpsf.add(ap.sum()*dA[id]);
	      avgpsf.push_back(aavgpsf.sum());

	      //std::cout << r2d(d[id]) << ' ' << ap.sum()*dA[id] << '\n';
	    }

	  std::vector<double>::iterator iavgpsf = 
	    std::upper_bound(avgpsf.begin(),avgpsf.end(),aavgpsf.sum()*0.68);
	  
	  std::cout 
	    << ' ' << r2d(d[iavgpsf-avgpsf.begin()]);
	}
      std::cout << '\n';
    }

  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<IRFs*>::iterator iirf=irfs.begin();iirf!=irfs.end();iirf++)
    delete *iirf;
}
