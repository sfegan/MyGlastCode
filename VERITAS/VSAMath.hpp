//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAMath.hpp
  Various math functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSAMATH_HPP
#define VSAMATH_HPP

#define VSA_ASSERT

#include<cmath>
#include<cfloat>

#ifdef VSA_ASSERT
#include<cassert>
#endif

namespace VERITAS
{
  namespace VSAMath
  {

    // ------------------------------------------------------------------------
    // Roots of Quadratic and Cubic
    // ------------------------------------------------------------------------

    unsigned realRoots2(const double a, const double b, const double c,
			double roots[2]);
    unsigned realRoots3(const double _a, const double _b, 
			const double _c, const double _d,
			double roots[3]);

    // ------------------------------------------------------------------------
    // Carlson forms of integrals
    // ------------------------------------------------------------------------

    // Carlson's elliptic integral of the first kind
    double ellipticRFCarlson(double x, double y, double z);

    inline double ellipticRFCarlsonZZero(const double x, const double y)
    {
      const double ERRTOL = 2.7*sqrt(DBL_EPSILON);
      double xm = sqrt(x);
      double ym = sqrt(y);
      while(fabs(xm-ym) > ERRTOL*fabs(xm))
	{
	  const double xm_next = 0.5*(xm+ym);
	  ym = sqrt(xm*ym);
	  xm = xm_next;
	}
      return M_PI/(xm+ym);
    }

    // Carlson's elliptic integral of the second kind
    double ellipticRDCarlson(double x, double y, double z);

    inline double ellipticRGCarlsonZZero(const double x, const double y)
    {
      const double ERRTOL = 2.7*sqrt(DBL_EPSILON);
      double xm = sqrt(x);
      double ym = sqrt(y);
      double fac = 0.5;
      double sum = (xm+ym)/2;
      sum *= sum;
      while(fabs(xm-ym) > ERRTOL*fabs(xm))
	{
	  const double xm_next = 0.5*(xm+ym);
	  ym = sqrt(xm*ym);
	  xm = xm_next;
	  const double dxy = xm-ym;
	  sum -= fac*dxy*dxy;
	  fac *= 2.0;
	}
      return 0.5*sum*M_PI/(xm+ym);
    }

    // ------------------------------------------------------------------------
    // Jacobi elliptic functions
    // ------------------------------------------------------------------------

    void ellipticJacobiFunctions(const double _u, const double k,
				 double &sn, double& cn, double& dn);

    // ------------------------------------------------------------------------
    // Elliptic integrals of the first kind
    // ------------------------------------------------------------------------


    // Legendre elliptic integral of the first kind
    inline double elliptic1Legendre(const double phi, const double k)
    {
      const double s    = sin(phi); 
      const double c    = cos(phi);      
      const double sskk = s*s*k*k;
      return s*ellipticRFCarlson(c*c,1.0-sskk,1.0);
    }

    // Complete elliptic integral of the first kind
    inline double elliptic1Complete(const double k)
    { 
      return ellipticRFCarlsonZZero(1.0-k*k,1.0);
    }

    // ------------------------------------------------------------------------
    // Elliptic integrals of the second kind
    // ------------------------------------------------------------------------

    // Legendre elliptic integral of the second kind
    inline double elliptic2Legendre(const double phi, const double k)
    {
      const double s    = sin(phi);
      const double c    = cos(phi);
      const double cc   = c*c;
      const double sskk = s*s*k*k;
      const double q    = 1.0-sskk;
      return s*(ellipticRFCarlson(cc,q,1.0) -
		sskk*ellipticRDCarlson(cc,q,1.0)/3.0);
    }

    // Complete elliptic integral of the second kind
    inline double elliptic2Complete(const double k)
    {
#if 0
      const double k2 = k*k;
      const double q  = 1.0-k2;
      return (ellipticRFCarlson(0.0,q,1.0)
	      - k2*ellipticRDCarlson(0.0,q,1.0)/3.0);
#endif
      return 2.0*ellipticRGCarlsonZZero(1.0-k*k,1.0);
    }

  }
}

#endif // VSAALGEBRA_HPP
