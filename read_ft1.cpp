//-*-mode:c++; mode:font-lock;-*-

#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "VERITAS/VSAAlgebra.hpp"

#include "FT1.hpp"
#include "FT2.hpp"
#include "Catalog.hpp"
#include "Util.hpp"
#include "MyMinuit.hpp"
#include "PSF.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

int main(int argc, char** argv)
{
  std::string progname(*argv);
  argc--, argv++;

  if(argc != 3)
    {
      std::cerr
	<< "Usage: " << progname << " ft1_file ft2_file catalog_file\n";
      exit(EXIT_FAILURE);
    }

  std::string event_file(*argv);
  argc--, argv++;

  std::string sc_file(*argv);
  argc--, argv++;

  std::string cat_file(*argv);
  argc--, argv++;

  std::vector<Catalog> catalog_vec;
  Catalog::loadFromFITS(catalog_vec, cat_file);

  Catalog::Filter cuts;
  cuts.glat_min = std::pair<bool,double>(true, d2r(18.5));
  cuts.ts_min = std::pair<bool,double>(true, 9);
  cuts.nassoc_min = std::pair<bool,unsigned>(true, 1);
  cuts.max_inter_assoc_dist_max = std::pair<bool,double>(true, d2r(0.1));
  cuts.min_neighbour_dist_min = std::pair<bool,double>(true, d2r(1.5));

  unsigned ncat = Catalog::filterCatalog(catalog_vec, cuts);
  std::cerr << "Catalog has: " << ncat << " entries\n";
  // Catalog::printCatalog(catalog_vec, std::cout, false);

  std::vector<double> cat_ra(ncat);
  std::vector<double> cat_dec(ncat);
  std::vector<double> cat_sra(ncat);
  std::vector<double> cat_cra(ncat);
  std::vector<double> cat_sdec(ncat);
  std::vector<double> cat_cdec(ncat);

#if 0
  for(unsigned icat=0;icat<ncat;icat++)
    {
      cat_ra[icat] = catalog_vec[icat].ra;
      cat_dec[icat] = catalog_vec[icat].dec;
    }
#else
  for(unsigned icat=0;icat<ncat;icat++)
    {
      cat_ra[icat] = catalog_vec[icat].assoc[0].ra;
      cat_dec[icat] = catalog_vec[icat].assoc[0].dec;
    }
#endif

  for(unsigned icat=0;icat<ncat;icat++)
    {
      cat_sra[icat] = std::sin(cat_ra[icat]);
      cat_cra[icat] = std::cos(cat_ra[icat]);
      cat_sdec[icat] = std::sin(cat_dec[icat]);
      cat_cdec[icat] = std::cos(cat_dec[icat]);
    }

  //std::vector<FT2> sc_vec;
  //  FT2::loadFromFITS(sc_vec, sc_file);

  std::vector<FT1> event_vec;
  unsigned nevent = FT1::loadFromFITS(event_vec, event_file);

  const double dcut = 1.5;
  std::vector<double> d_vec;
  for(unsigned ievent=0;ievent<nevent;ievent++)
    {
      FT1& event(event_vec[ievent]);
      double ra = event.ra;
      double dec = event.dec;

      double sra = std::sin(ra);
      double cra = std::cos(ra);
      double sdec = std::sin(dec);
      double cdec = std::cos(dec);

      double dmin = 0;
      unsigned dmin_icat = 0;
      for(unsigned icat=0;icat<ncat;icat++)
	{
	  double d = 
	    approx_sphere_dist_nt(sra,cra,sdec,cdec,
				  cat_sra[icat],cat_cra[icat],
				  cat_sdec[icat],cat_cdec[icat]);
	  if(icat==0 || d<dmin)dmin=d,dmin_icat=icat;
	}

      dmin = sphere_dist(ra,dec,cat_ra[dmin_icat],cat_dec[dmin_icat]);
      dmin = r2d(dmin);
      if(dmin<dcut)d_vec.push_back(dmin);
    }

  double r68_g;
  if(1)
    {
      Like_GaussianPSFWithBkg psf(d_vec, dcut);
      MyMinuit minuit(&psf, 1);
      minuit.minimize();
      
      unsigned np = psf.numParam();
      for(unsigned ip=0;ip<np;ip++)
	std::cout << minuit.pOpt(ip) << ' ' << minuit.dpOpt(ip) << '\n';
      
      r68_g = minuit.pOpt(0)*std::sqrt(-2*log(1-0.68));
    }
  
  double r68_b;
  if(1)
    {
      Like_BurdettianPSFWithBkg psf(d_vec, dcut);
      MyMinuit minuit(&psf, 1);
      minuit.minimize();
      
      unsigned np = psf.numParam();
      for(unsigned ip=0;ip<np;ip++)
	std::cout << minuit.pOpt(ip) << ' ' << minuit.dpOpt(ip) << '\n';
      
      double gamma = minuit.pOpt(1);
      r68_b = 
	minuit.pOpt(0)*
	std::sqrt(2*gamma*(std::pow(1.0-0.68,-1.0/(gamma-1.0))-1.0));
    }

  std::cout << "r86_g = " << r68_g << '\n';
  std::cout << "r86_b = " << r68_b << '\n';
}
