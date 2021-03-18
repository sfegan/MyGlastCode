#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "RegularSpline.hpp"
#include "SourceLibrary.hpp"

int main(int argc, char** argv)
{
  std::string progname(*argv);
  argc--,argv++;

  if(argc != 2)
    {
      std::cerr << "Usage: " << progname << " filename.xml source_name\n";
      std::exit(EXIT_FAILURE);
    }

  std::string filename(*argv);
  argc--,argv++;

  std::string srcname(*argv);
  argc--,argv++;

  SourceLibrary slib;
  slib.readFromXML(filename);

  const SourceLibrary::Source* s(slib.source(srcname));
  if(s && s->spectrum().type() == "SplineFunction")
    {
      unsigned npoints = 0;
      if(s->spectrum().attributeAs("npoints",npoints) == false)
	{
	  std::cerr << "Could not get number of spline points!\n";
	  std::exit(EXIT_FAILURE);
	}
      
      const SourceLibrary::Parameter* p = 0;

      double xmin = 100;
      p = s->spectrum().parameter("LowerLimit");
      if(p)xmin = p->value();
      xmin = std::log10(xmin);

      double xmax = 100000;
      p = s->spectrum().parameter("UpperLimit");
      if(p)xmax = p->value();
      xmax = std::log10(xmax);

      double gamma1 = -2.0;
      double dgamma1 = 0;
      p = s->spectrum().parameter("Index1");
      if(p)gamma1 = p->value(), dgamma1 = p->error();

      double gammaN = -2.0;
      double dgammaN = 0;
      p = s->spectrum().parameter("IndexN");
      if(p)gammaN = p->value(), dgammaN = p->error();
      
      std::vector<double> yi(npoints);
      std::vector<double> dyi(npoints);
      for(unsigned ip=0;ip<npoints;ip++)
	{
	  char buffer[80];
	  sprintf(buffer,"Log10Flx%d",ip);
	  p = s->spectrum().parameter(buffer);
	  if(p)yi[ip] = p->value(), dyi[ip] = p->error();
	}

      RegularSpline spline(xmin,xmax,yi,
			   RegularSpline::BCS_CLAMPED,gamma1,
			   RegularSpline::BCS_CLAMPED,gammaN,
			   true);
      
      int ni = (npoints-1)*100;
      for(int i=0;i<=ni;i++)
	{
	  double x = double(i)/double(ni)*(xmax-xmin) + xmin;
	  double f = spline.f(x);
	  double df2 = 0;
	  double dfc = 0;

	  dfc = spline.dfdSi(x,RegularSpline::S_1)*dgamma1;
	  df2 += dfc*dfc;

	  dfc = spline.dfdSi(x,RegularSpline::S_N)*dgammaN;
	  df2 += dfc*dfc;

	  for(unsigned ip=0;ip<npoints;ip++)
	    {
	      dfc = spline.dfdYi(x,ip)*dyi[ip];
	      df2 += dfc*dfc;
	    }
	  
	  double df = std::sqrt(df2);
	  
	  std::cout << x << ' ' << std::pow(10,x) << ' ' 
		    << f << ' ' << df << ' '
		    << std::pow(10,f) << ' ' << df*std::pow(10,f) << ' '
		    << std::pow(10,f+df) << ' ' << std::pow(10,f-df) << ' '
		    << ((i%100)==0?1:0) << '\n';
	}
#if 0
#endif
    }
  else if(s)
    {
      std::cerr << "Spectrum for source \"" << srcname 
		<< "\" is not SplineFunction (\"" << s->spectrum().type()
		<< "\")\n";
      std::exit(EXIT_FAILURE);
    }
  else
    {
      std::cerr << "Could not get spectrum for source \"" << srcname << "\"\n";
      std::exit(EXIT_FAILURE);
    }

}
