//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAMath.cpp
  Various math functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#include<cmath>
#include<cfloat>
#include<limits>
#include<algorithm>

#include<VSAMath.hpp>

using namespace VERITAS;

static inline double sgn(double x)
{
  if(x>=0)return 1.0;
  else return -1.0;
}

unsigned VSAMath::realRoots2(const double a, const double b, const double c,
			     double roots[2])
{
  // Reference: Numerical Recipies

  // a - coefficient of x^2
  // b - coefficient of x^1
  // c - coefficient of x^0

  // formula x = ( -b +/- sqrt(b^2 - 4 a c) ) / 2 a

#ifdef VSA_ASSERT
  assert(a!=0);
#endif

  const double R = b*b - 4.0*a*c;
  
  if(R<0)
    {
      roots[0] = roots[1] = FP_NAN;
      return 0;
    }
  else if(R==0)
    {
      roots[0] = roots[1] = -b/2.0/a;
      return 1;
    }
  else
    {
      double q = -0.5 * ( b + sgn(b)*sqrt(R) );
      roots[0] = q/a;
      roots[1] = c/q;
      if(roots[0]>roots[1])std::swap(roots[0],roots[1]);
      return 2;
    }
  
#ifdef VSA_ASSERT
  assert(0); // never reach here
#endif
}

unsigned VSAMath::realRoots3(const double _a, const double _b, 
			     const double _c, const double _d,
			     double roots[3])
{
  // Reference: Numerical Recipies

  // _a - coefficient of x^3
  // _b - coefficient of x^2
  // _c - coefficient of x^1
  // _d - coefficient of x^0

#ifdef VSA_ASSERT
  assert(_a!=0);
#endif

  // NR uses different variables -- translate for clarity
  const double a = _b/_a;
  const double b = _c/_a;
  const double c = _d/_a;

  const double a2 = a*a;
  const double a3 = a2*a;

  const double Q = (a2 - 3.0*b)/9.0;
  const double R = (2.0*a3 - 9.0*a*b + 27.0*c)/54.0;

  const double R2 = R*R;
  const double Q3 = Q*Q*Q;

  double decider = R2-Q3;
  if(fabs(R2+Q3)>0)decider /= fabs(R2+Q3);

  if(decider < -DBL_EPSILON)
    {
      // Three real roots
      const double theta = acos(R/sqrt(Q3));
      const double factor0 = -2.0*sqrt(Q);
      const double factor1 = a/3.0;
      roots[0] = factor0 * cos(theta/3.0) - factor1;
      roots[1] = factor0 * cos((theta+2.0*M_PI)/3.0) - factor1;
      roots[2] = factor0 * cos((theta-2.0*M_PI)/3.0) - factor1;
      if(roots[0]>roots[1])std::swap(roots[0],roots[1]);
      if(roots[1]>roots[2])
	{
	  std::swap(roots[1],roots[2]);
	  if(roots[0]>roots[1])std::swap(roots[0],roots[1]);
	}
      return 3;
    }
  else if(decider < DBL_EPSILON)
    {
      if(R==0)
	{
	  // One real root
	  roots[0] = roots[1] = roots[2] = -a/3.0;
	  return 1;
	}
      else
	{
	  const double A = -sgn(R)*pow(fabs(R),1.0/3.0);
	  // Two real roots
	  roots[0] = 2.0*A-a/3.0;
	  roots[1] = roots[2] = -A-a/3.0;
	  if(roots[0]>roots[2])std::swap(roots[0],roots[2]);
	  return 2;
	}
    }
  else
    {
      const double A = -sgn(R)*pow(fabs(R) + sqrt(R2-Q3),1.0/3.0);
      const double B = (A==0)?0:Q/A;
      roots[0] = (A+B)-a/3.0;

      if(fabs(A-B)<DBL_EPSILON)
	{
	  //.I think this can never happen 

	  // A-B=0 => R^2-Q^3=0  /OR/  R^2-Q^3=|2R|
	  // The first case is manifestly not true in this case (decider>0)
	  // Is the second possible -- have to multiply it out some time

	  roots[1] = roots[2] = -0.5*(A+B)-a/3.0;
	  if(roots[0]>roots[2])std::swap(roots[0],roots[2]);
	  return 2;
	}
      else
	{
	  roots[1] = roots[2] = FP_NAN;
	  return 1;
	}
    }

#ifdef VSA_ASSERT
  assert(0); // never reach here
#endif
}

// ----------------------------------------------------------------------------
// Carlson forms of integrals
// ----------------------------------------------------------------------------

double VSAMath::ellipticRFCarlson(double x, double y, double z)
{
  static const double ERRTOL = 0.0025;
  static const double TINY   = // 1.5e-38
    1.1*std::numeric_limits<double>::min()*5.0; 
  static const double BIG    = // 3.0e37
    0.9*std::numeric_limits<double>::max()*0.2;
  static const double THIRD  = 1.0/3.0;
  static const double C1     = 1.0/24.0;
  static const double C2     = 0.1;
  static const double C3     = 3.0/44.0;
  static const double C4     = 1.0/14.0;

  assert((x>=0)&&(y>=0)&&(z>=0)&&(x<=BIG)&&(y<=BIG)&&(z<=BIG)
	 &&(std::min(std::min(x+y,x+z),y+z)>=TINY));

  double xt = x;
  double yt = y;
  double zt = z;

  double delx;
  double dely;
  double delz;
  double ave;

  do {
    const double sqrtx = sqrt(xt);
    const double sqrty = sqrt(yt);
    const double sqrtz = sqrt(zt);
    const double alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    xt = 0.25*(xt+alamb);
    yt = 0.25*(yt+alamb);
    zt = 0.25*(zt+alamb);
    ave = THIRD*(xt+yt+zt);
    delx = (ave-xt)/ave;
    dely = (ave-yt)/ave;
    delz = (ave-zt)/ave;
  } while((fabs(delx)>ERRTOL)||(fabs(dely)>ERRTOL)||(fabs(delz)>ERRTOL));

  const double e2 = delx*dely - delz*delz;
  const double e3 = delx*dely*delz;
  return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}

// Carlson's elliptic integral of the second kind
double VSAMath::ellipticRDCarlson(double x, double y, double z)
{
  static const double ERRTOL = 0.0015;
  static const double TINY   = //1e-25;
    1.1*pow(std::numeric_limits<double>::max(),-2.0/3.0)*2.0;
  static const double BIG    = //4.5e21;
    0.9*pow(std::numeric_limits<double>::min(),-2.0/3.0)*0.1*ERRTOL;
  static const double C1     = 3.0/14.0;
  static const double C2     = 1.0/6.0;
  static const double C3     = 9.0/22.0;
  static const double C4     = 3.0/26.0;
  static const double C5     = 0.25*C3;
  static const double C6     = 1.5*C4;

  assert((x>=0)&&(y>=0)&&(x+y>=TINY)&&(z>=TINY)&&(x<=BIG)&&(y<=BIG)&&(z<=BIG));
  
  double xt = x;
  double yt = y;
  double zt = z;
  double sum = 0.0;
  double fac = 1.0;

  double delx;
  double dely;
  double delz;
  double ave;

  do {
    const double sqrtx = sqrt(xt);
    const double sqrty = sqrt(yt);
    const double sqrtz = sqrt(zt);
    const double alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    sum += fac/(sqrtz*(zt+alamb));
    fac  = 0.25*fac;
    xt   = 0.25*(xt+alamb);
    yt   = 0.25*(yt+alamb);
    zt   = 0.25*(zt+alamb);
    ave  = 0.2*(xt+yt+3.0*zt);
    delx = (ave-xt)/ave;
    dely = (ave-yt)/ave;
    delz = (ave-zt)/ave;
  } while((fabs(delx)>ERRTOL)||(fabs(dely)>ERRTOL)||(fabs(delz)>ERRTOL));
  
  const double ea = delx*dely;
  const double eb = delz*delz;
  const double ec = ea-eb;
  const double ed = ea-6.0*eb;
  const double ee = ed+ec+ec;

  return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		      + delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}

// ----------------------------------------------------------------------------
// Jacobi elliptic functions
// ----------------------------------------------------------------------------
#include<iostream>
void VSAMath::ellipticJacobiFunctions(const double _u, const double k,
				      double &sn, double& cn, double& dn)
{
  static const double CA = 2.7*sqrt(DBL_EPSILON);

  double emc = 1.0-k*k;

  if(emc) 
    {
      bool bo = (emc < 0.0);
      double u = _u;
      double d = 0;
      if(bo)
	{
	  d = 1.0-emc;
	  emc = -d;
	  d = sqrt(d);
	  u *= d;
	}


      dn = 1.0;

      double a = 1.0;
      double c;
      double em[13];
      double en[13];
      int l;
      for(int i=0;i<13;i++)
	{
	  l = i;
	  em[i] = a;
	  emc = sqrt(emc);
	  en[i] = emc;
	  c = 0.5*(a+emc);
	  if(fabs(a-emc) <= CA*a)break;
	  emc *= a;
	  a = c;
	}
      u *= c;

      sn = sin(u);
      cn = cos(u);

      if(sn) 
	{
	  a = cn/sn;
	  c *= a;
	  for(int i=l;i>=0;i--)
	    {
	      double b = em[i];
	      a *= c;
	      c *= dn;
	      dn = (en[i]+a)/(b+a);
	      a = c/b;
	    }
	  a = 1.0/sqrt(c*c+1.0);
	  sn=(sn >= 0.0 ? a : -a);
	  cn=c*sn;
	}
      if (bo)
	{
	  a = dn;
	  dn = cn;
	  cn = a;
	  sn /= d;
	}
    } 
  else 
    {
      cn = 1.0/cosh(_u);
      dn = cn;
      sn = tanh(_u);
    }
}


#ifdef TESTMAIN
#include<iostream>
int main(int argc, char** argv)
{
  for(double x=0;x<1.0;x+=0.01)
    std::cout << x << ' ' 
	      << VSAMath::elliptic1Complete(x) << ' ' 
	      << VSAMath::elliptic2Complete(x) << '\n'; 
}
#endif
