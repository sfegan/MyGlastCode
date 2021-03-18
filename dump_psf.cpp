//-*-mode:c++; mode:font-lock;-*-

// dump_psf.cpp - Stephen Fegan 
//              - sfegan@llr.in2p3.fr
//              - 2010-07-28
//
// Dump the PSF vs ENERGY from the IRFs
//
// $Id: lomblike.cpp 1995 2010-07-26 14:33:48Z sfegan $

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

  double theta = 0;
  options.findWithValue("theta", theta,
			"Set desired theta angle.");
  theta = d2r(theta);

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

  double costheta=cos(theta);
  for(double log10e=1;log10e<6.001;log10e+=0.05)
    {
      std::cout << log10e << ' ' << std::pow(10.0,log10e);
      for(unsigned iirf=0;iirf<irfs.size();iirf++)
	{
	  double dd = 0.001;
	  Accumulator axx;
	  Accumulator a;
	  for(double d=0.5*dd;d<10;d+=dd)
	    {
	      double dA = 
		d2r(d2r(M_PI*((d+dd*0.5)*(d+dd*0.5)-(d-dd*0.5)*(d-dd*0.5))));
	      const double p = irfs[iirf]->psf()->value(d2r(d),log10e,costheta);
	      //std::cout << d << ' ' << p << '\n';
	      a.add(p*dA);
	      axx.add(p*d*d*dA);
	    }
	  std::cout 
	    << ' ' << sqrt(axx.sum()/a.sum()) << ' ' << a.sum() 
	    << ' ' << irfs[iirf]->psf()->integral(d2r(10),log10e,costheta)
	    << ' ' << r2d(irfs[iirf]->psf()->rP(0.68,log10e,costheta));
	}
      std::cout << '\n';
    }

  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<IRFs*>::iterator iirf=irfs.begin();iirf!=irfs.end();iirf++)
    delete *iirf;
}
