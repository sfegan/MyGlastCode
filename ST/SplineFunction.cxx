/** 
 * @file SplineFunction.cxx
 * @brief Implementation for the SplineFunction Function class
 * @author S. Fegan
 *
 * $Header: $
 */

#include<vector>
#include<string>
#include<limits>
#include<cmath>
#include<exception>

#include "optimizers/dArg.h"
#include "optimizers/Exception.h"
#include "optimizers/ParameterNotFound.h"
#include "Likelihood/SplineFunction.h"

namespace Likelihood {

// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================

/* 
  
   RegularSpline.hpp - Stephen Fegan 
                     - sfegan@llr.in2p3.fr
                     - 19 February 2010

   Class to evaluate regular spline and derivatives with respect to
   parameters. For information on cubic splines and how to solve
   tridiagonal matrices consult:

   http://mathworld.wolfram.com/CubicSpline.html
   http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
   Numerical recipes (3rd edition, section 3.3)

   $Id: $

*/

#include<vector>
#include<cmath>

class RegularSpline
{
public:
  enum BndryCond { BCS_NATURAL, BCS_NOT_A_KNOT, BCS_CLAMPED };
  enum Slope { S_1, S_N };

  RegularSpline(double xmin, double xmax, const std::vector<double>& yi,
	        BndryCond bcs1 = BCS_NOT_A_KNOT, double slope1 = 0,
		BndryCond bcsN = BCS_NOT_A_KNOT, double slopeN = 0,
		bool enable_derivatives = false);

  double f(double x) const { double y; f(x,y); return y; }
  double dfdx(double x) const { double y; dfdx(x,y); return y; }
  double d2fdx2(double x) const { double y; d2fdx2(x,y); return y; }
  double d3fdx3(double x) const { double y; d3fdx3(x,y); return y; }
  double dfdYi(double x, unsigned ip) const
  { double dydyi; dfdYi(x,ip,dydyi); return dydyi; }
  double dfdSi(double x, Slope is) const
  { double dydsi; dfdSi(x,is,dydsi); return dydsi; }

  bool f(double x, double& y) const;
  bool dfdx(double x, double& dydx) const;
  bool d2fdx2(double x, double& d2ydx2) const;
  bool d3fdx3(double x, double& d3ydx3) const;
  bool dfdYi(double x, unsigned ip, double& dydyi) const;
  bool dfdSi(double x, Slope is, double& dydsi) const;

private:
  struct SplinePiece
  {
    double              x1;         // x at the startpoint of the spline
    double              x2;         // x at the endpoint of the spline
    double              y1;         // y at the start of the spline
    double              y2;         // y at the end of the spline
    double              y1pp6;      // (d2y/dx^2)/6 at the start of the spline
    double              y2pp6;      // (d2y/dx^2)/6 at the end of spline
    std::vector<double> dy1pp6_dyi; // d[(d2y/dx^2)/6]/dYi at start of spline
    std::vector<double> dy2pp6_dyi; // d[(d2y/dx^2)/6]/dYi at end of spline
    double              dy1pp6_ds1; // d[(d2y/dx^2)/6]/ds1 at start of spline
    double              dy2pp6_ds1; // d[(d2y/dx^2)/6]/ds1 at end of spline
    double              dy1pp6_dsN; // d[(d2y/dx^2)/6]/dsN at start of spline
    double              dy2pp6_dsN; // d[(d2y/dx^2)/6]/dsN at end of spline
  };

  unsigned SAB(double x, const SplinePiece*& spc, double& A, double& B) const
  {
    unsigned interval_no;
    if(x>=m_xNthBinLHS)interval_no=m_nspline-1;
    else if(x<m_x1stBinRHS)interval_no=0;
    else interval_no=floorl((x-m_xmin)/m_dx);

    spc = &m_spline[interval_no];
    A = (spc->x2 - x)/m_dx;
    B = 1-A;
    return interval_no;
  }

  bool                       m_enable_derivatives;
  double                     m_xmin;
  double                     m_xmax;
  double                     m_dx;
  double                     m_x1stBinRHS;
  double                     m_xNthBinLHS;
  double                     m_nspline;
  std::vector<SplinePiece>   m_spline;
};

// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================

/* 
  
   RegularSpline.cpp - Stephen Fegan 
                     - sfegan@llr.in2p3.fr
                     - February 2010

   Class to evaluate regular spline and derivatives with respect to
   parameters. For information on cubic splines and how to solve
   tridiagonal matrices consult:

   http://mathworld.wolfram.com/CubicSpline.html
   http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorith
   Numerical recipes (3rd edition, section 3.3)

   $Id: $

*/

#include<string>

RegularSpline::
RegularSpline(double xmin, double xmax, const std::vector<double>& yi,
	      BndryCond bcs1, double slope1,
	      BndryCond bcsN, double slopeN,
	      bool enable_derivatives):
  m_enable_derivatives(enable_derivatives),
  m_xmin(xmin), m_xmax(xmax), m_dx((xmax-xmin)/double(yi.size()-1)), 
  m_x1stBinRHS(), m_xNthBinLHS(),
  m_nspline(yi.size()-1), m_spline(m_nspline)
{
  unsigned np = yi.size();
  
  // First solve the tridiagonal matrix equation for the slopes at
  // each of the points.. do this in a manner so that we can later
  // calculate the value of the spline and the derivative with respect
  // to the input points and slopes

  std::vector<double> a(np, 1); // Sub diagnoal elements   [1,np-1]
  std::vector<double> b(np, 4); // Diagonal elements       [0,np-1]
  std::vector<double> c(np, 1); // Super diagonal elements [0,np-2]
  std::vector<double> d(np, 0); // RHS of equation         [0,np-1]
  std::vector<double> m(np, 0); // Slopes we are seeking   [0,np-1]

  std::vector<std::vector<double> > dp;
  std::vector<double> ds1;
  std::vector<double> dsN;
  std::vector<std::vector<double> > mp;
  std::vector<double> ms1;
  std::vector<double> msN;

  if(m_enable_derivatives)
    {
      dp.resize(np, std::vector<double>(np,0));
      ds1.resize(np, 0);
      dsN.resize(np, 0);
      mp.resize(np, std::vector<double>(np,0));
      ms1.resize(np, 0);
      msN.resize(np, 0);
    }

  // We cannot use the default NOT_A_KNOT algorithm with only 2 or 3
  // points (i.e. 1 or 2 spline segments), since there is no "knot" in
  // the first case, giving no constraints at all, and only one "knot"
  // in the second case giving a single constraint. 

  // In the 2-point case, if either of the BCS is set to NOT_A_KNOT we
  // zero out the cubic term. If both are set to zero then we also set
  // the quadratic term to zero (by changing bcsN to BCS_NATURAL). 
  // This has the pleasent effect of making the spline into a line
  // segment if either is NOT_A_KNOT and the other is either NATURAL
  // or NOT_A_KNOT. In the case of NOT_A_KNOT combined with CLAMPED 
  // (preassigned slope) we get a quadratic. 

  // In the 3-point case, if either of the BCS is set to NOT_A_KNOT we
  // set the quadratic term to be constant across the "knot" as
  // expected. The other constraint then depends on what the other
  // BCS is. For NATURAL and CLAMPED the other constraint is obvious.
  // If the other is also NOT_A_KNOT we set the cubic term to zero to
  // get a quadratic that fits the three points.
  
  switch(bcs1)
    {
    case BCS_NATURAL:
      b[0] = 2;
      d[0] = 3*(yi[1]-yi[0]);
      if(m_enable_derivatives)
	dp[0][0] = -3,
	  dp[0][1] = 3;
      break;
    case BCS_CLAMPED:
      b[0] = 1;
      c[0] = 0;
      d[0] = slope1*m_dx;
      if(m_enable_derivatives)
	ds1[0] = m_dx;
      break;
    case BCS_NOT_A_KNOT:
      if((np==2)||((np==3)&&(bcsN==BCS_NOT_A_KNOT)))
	{
	  // Set cubic term to zero
	  b[0] = 1;
	  d[0] = 2*(yi[1]-yi[0]);
	  if(m_enable_derivatives)
	    dp[0][0] = -2,
	      dp[0][1] = 2;
	  if((np==2)&&(bcsN==BCS_NOT_A_KNOT))bcsN = BCS_NATURAL;
	}
      else
	{
	  b[0] = 2;
	  c[0] = 4;
	  d[0] = yi[2]+4*yi[1]-5*yi[0];
	  if(m_enable_derivatives)
	    dp[0][0] = -5,
	      dp[0][1] = 4,
	      dp[0][2] = 1;
	}
      break;
    }

  switch(bcsN)
    {
    case BCS_NATURAL:
      b[np-1] = 2;
      d[np-1] = 3.0*(yi[np-1]-yi[np-2]);
      if(m_enable_derivatives)
	dp[np-1][np-2] = -3,
	  dp[np-1][np-1] = 3;
      break;
    case BCS_CLAMPED:
      b[np-1] = 1;
      a[np-1] = 0;
      d[np-1] = slopeN*m_dx;
      if(m_enable_derivatives)
	dsN[np-1] = m_dx;
      break;
    case BCS_NOT_A_KNOT:
      if(np==2)
	{
	  // Set cubic term = 0
	  b[1] = 1;
	  d[1] = 2*(yi[1]-yi[0]);
	  if(m_enable_derivatives)
	    dp[1][0] = -2,
	      dp[1][1] = 2;
	}
      else
	{
	  b[np-1] = 2;
	  a[np-1] = 4;
	  d[np-1] = -yi[np-3]-4.0*yi[np-2]+5*yi[np-1];
	  if(m_enable_derivatives)
	    dp[np-1][np-3] = -1, 
	      dp[np-1][np-2] = -4, 
	      dp[np-1][np-1] = 5;
	}
      break;
    }
  
  for(unsigned ip=1; ip<np-1; ip++)
    d[ip] = 3.0*(yi[ip+1]-yi[ip-1]);

  if(m_enable_derivatives)
    for(unsigned ip=1; ip<np-1; ip++)
      dp[ip][ip-1] = -3,
	dp[ip][ip+1] = 3;
  
#if 0
  for(unsigned ip=0; ip<np; ip++)
    std::cout << a[ip] << ' ';
  std::cout << '\n';
  for(unsigned ip=0; ip<np; ip++)
    std::cout << b[ip] << ' ';
  std::cout << '\n';
  for(unsigned ip=0; ip<np; ip++)
    std::cout << c[ip] << ' ';
  std::cout << '\n';
  for(unsigned ip=0; ip<np; ip++)
    std::cout << d[ip] << ' ';
  std::cout << '\n';
#endif

  c[0] /= b[0];
  d[0] /= b[0];
  if(m_enable_derivatives)
    {
      for(unsigned jp=0;jp<np;jp++)
	dp[0][jp] /= b[0];
      ds1[0] /=  b[0];
      dsN[0] /=  b[0];
    }
  for(unsigned ip = 1; ip<np; ip++)
    {
      double id = 1.0/(b[ip] - c[ip-1]*a[ip]);
      c[ip] *= id;
      d[ip] = (d[ip] - d[ip-1] * a[ip])*id;
      if(m_enable_derivatives)
	{
	  for(unsigned jp=0;jp<np;jp++)
	    dp[ip][jp] = (dp[ip][jp] - dp[ip-1][jp] * a[ip])*id;
	  ds1[ip] = (ds1[ip] - ds1[ip-1] * a[ip])*id;
	  dsN[ip] = (dsN[ip] - dsN[ip-1] * a[ip])*id;
	}
    }
 
  m[np-1] = d[np-1];
  for (int ip=np-2; ip>=0; ip--)
    m[ip] = d[ip]-c[ip]*m[ip+1];

  if(m_enable_derivatives)
    {
      for(unsigned jp=0;jp<np;jp++)
	mp[np-1][jp] = dp[np-1][jp];
      ms1[np-1] = ds1[np-1];
      msN[np-1] = dsN[np-1];

      for(int ip=np-2; ip>=0; ip--)
	{
	  for(unsigned jp=0;jp<np;jp++)
	    mp[ip][jp] = dp[ip][jp]-c[ip]*mp[ip+1][jp];
	  ms1[ip] = ds1[ip]-c[ip]*ms1[ip+1];
	  msN[ip] = dsN[ip]-c[ip]*msN[ip+1];
	}
    }

#if 0
  for(unsigned ip=0; ip<np; ip++)
    std::cout << m[ip] << ' ';
  std::cout << '\n';
#endif

  // Now prepare the spline structure
  
  for(unsigned ispline=0; ispline<m_nspline; ispline++)
    {
      double i1 = double(ispline)/double(m_nspline);
      double i2 = double(ispline+1)/double(m_nspline);
      double x1 = m_xmin*(1.0-i1) + m_xmax*i1;
      double x2 = m_xmin*(1.0-i2) + m_xmax*i2;
      double y1 = yi[ispline];                    // y-value at x1
      double y2 = yi[ispline+1];                  // y-value at x2
      double m1 = m[ispline];                     // first derivative at x1
      double m2 = m[ispline+1];                   // first derivative at x2

      double dx = x2-x1;
      double dy = y2-y1;

      double y1pp = 2*( 3*dy - 2*m1 - m2);        // second derivative at x1
      double y2pp = 2*(-3*dy + m1 + 2*m2);        // second derivative at x2

      if(fabs(dx-m_dx)/m_dx >1e-6)
	throw std::string("Spline not regular");
      
      m_spline[ispline].x1    = x1;
      m_spline[ispline].x2    = x2;
      m_spline[ispline].y1    = y1;
      m_spline[ispline].y2    = y2;
      m_spline[ispline].y1pp6 = y1pp/6;
      m_spline[ispline].y2pp6 = y2pp/6;

      if(m_enable_derivatives)
	{
	  m_spline[ispline].dy1pp6_dyi.resize(np,0);
	  m_spline[ispline].dy2pp6_dyi.resize(np,0);

	  m_spline[ispline].dy1pp6_dyi[ispline]   += 2*(-3.0)/6;
	  m_spline[ispline].dy2pp6_dyi[ispline]   += 2*( 3.0)/6;
	  m_spline[ispline].dy1pp6_dyi[ispline+1] += 2*( 3.0)/6;
	  m_spline[ispline].dy2pp6_dyi[ispline+1] += 2*(-3.0)/6;

	  for(unsigned ip=0;ip<np;ip++)
	    {
	      double dm1dyi = mp[ispline][ip];
	      double dm2dyi = mp[ispline+1][ip];

	      m_spline[ispline].dy1pp6_dyi[ip]    += 2*(-2*dm1dyi-dm2dyi)/6;
	      m_spline[ispline].dy2pp6_dyi[ip]    += 2*( dm1dyi+2*dm2dyi)/6;
	    }
	  
	  double dm1ds1 = ms1[ispline];
	  double dm1dsN = msN[ispline];
	  double dm2ds1 = ms1[ispline+1];
	  double dm2dsN = msN[ispline+1];
	  m_spline[ispline].dy1pp6_ds1             = 2*(-2*dm1ds1-dm2ds1)/6;
	  m_spline[ispline].dy1pp6_dsN             = 2*(-2*dm1dsN-dm2dsN)/6;
	  m_spline[ispline].dy2pp6_ds1             = 2*( dm1ds1+2*dm2ds1)/6;
	  m_spline[ispline].dy2pp6_dsN             = 2*( dm1dsN+2*dm2dsN)/6;
	}
    }

  m_x1stBinRHS     = m_spline[0].x2;
  m_xNthBinLHS     = m_spline[m_nspline-1].x1;
}

bool RegularSpline::f(double x, double& y) const
{
  const SplinePiece* spc;
  double A,B;
  SAB(x, spc, A, B);
      
  y = A*(spc->y1+(A*A-1)*spc->y1pp6) + B*(spc->y2+(B*B-1)*spc->y2pp6);

  return true;
}

bool RegularSpline::dfdx(double x, double& dydx) const
{
  const SplinePiece* spc;
  double A,B;
  SAB(x, spc, A, B);
      
  dydx = (-(spc->y1+(3*A*A-1)*spc->y1pp6)+(spc->y2+(3*B*B-1)*spc->y2pp6))/m_dx;

  return true;
}

bool RegularSpline::d2fdx2(double x, double& d2ydx2) const
{
  const SplinePiece* spc;
  double A,B;
  SAB(x, spc, A, B);

  d2ydx2 = 6*(A*spc->y1pp6 + B*spc->y2pp6)/(m_dx*m_dx);
  
  return true;
}

bool RegularSpline::d3fdx3(double x, double& d3ydx3) const
{
  const SplinePiece* spc;
  double A,B;
  SAB(x, spc, A, B);

  d3ydx3 = 6*(-spc->y1pp6 + spc->y2pp6)/(m_dx*m_dx*m_dx);
  
  return true;
}

bool RegularSpline::dfdYi(double x, unsigned ip, double& dydyi) const
{
  if((!m_enable_derivatives)||(ip>m_nspline))return false;

  const SplinePiece* spc;
  double A,B;
  unsigned interval_no = SAB(x, spc, A, B);

  dydyi = A*(A*A-1)*spc->dy1pp6_dyi[ip] + B*(B*B-1)*spc->dy2pp6_dyi[ip];
  if(ip == interval_no)dydyi += A;
  else if(ip == interval_no+1)dydyi += B;
  
  return true;
}

bool RegularSpline::dfdSi(double x, Slope is, double& dydsi) const
{
  if(!m_enable_derivatives)return false;
  const SplinePiece* spc;
  double A,B;
  SAB(x, spc, A, B);

  double dy1pp6_ds;
  double dy2pp6_ds;
  if(is==S_1)dy1pp6_ds = spc->dy1pp6_ds1, dy2pp6_ds = spc->dy2pp6_ds1;
  else       dy1pp6_ds = spc->dy1pp6_dsN, dy2pp6_ds = spc->dy2pp6_dsN;
  
  dydsi = A*(A*A-1)*dy1pp6_ds + B*(B*B-1)*dy2pp6_ds;
  return true;
}

// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================

SplineFunction::
SplineFunction(double LowerLimit, double UpperLimit):
  Function(),
  m_initialized(false), m_npoints(), m_cpp(CPP_EDGE),
  m_cValues(), m_spline()
{
  typedef optimizers::Parameter Parameter;
  setMaxNumParams(3);

  addParam(Parameter("LowerLimit",    LowerLimit, 10,   1e6, false));
  addParam(Parameter("UpperLimit",    UpperLimit, 10,   1e6, false));
  addParam(Parameter("Normalization", 1.0,        1e-6, 1e6, false));
  
  setParamAlwaysFixed("LowerLimit");
  setParamAlwaysFixed("UpperLimit");

  // Set FuncType and ArgType for use with CompositeFunction hierarchy.
  m_funcType    = Addend;
  m_argType     = "dArg";
  m_genericName = "SplineFunction";
  m_normParName = "Normalization";  
}

double SplineFunction::value(optimizers::Arg & xarg) const
{
  updateCache();
  double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
  double log10X = std::log10(x);
  return m_parameter[Normalization].getTrueValue()*
    std::pow(10.0, m_spline->f(log10X));
}

double SplineFunction::derivByParam(optimizers::Arg & xarg, 
				    const std::string & paramName) const
{
  updateCache();
  double x = dynamic_cast<optimizers::dArg &>(xarg).getValue();
  double log10X = std::log10(x);

  int iparam(-1);
  for (unsigned int i = 0; i < m_parameter.size(); i++) {
    if (paramName == m_parameter[i].getName()) {
      iparam = i;
      break;
    }
  }

  if (iparam == -1) {
    throw optimizers::ParameterNotFound(paramName, getName(),
					"SplineFunction::derivByParam");
  }
  
  double dlog10FluxdP = 0;
  switch (iparam) 
    {
    case LowerLimit:
    case UpperLimit:
      throw optimizers::Exception("SplineFunction::derivByParam: attempt to "
				  "take derivative wrt a fixed parameter.");
      break;
    case Normalization:
      return m_parameter[Normalization].getScale()*
	std::pow(10.0, m_spline->f(log10X));
      break;
    case Index1:
      dlog10FluxdP = m_spline->dfdSi(log10X, RegularSpline::S_1);
      break;
    case IndexN:
      dlog10FluxdP = m_spline->dfdSi(log10X, RegularSpline::S_N);
      break;
    default:
      {
	unsigned ip = iparam - unsigned(Log10Flux0);
	dlog10FluxdP = m_spline->dfdYi(log10X, ip);
      }
      break;
    }
  
  return m_parameter[Normalization].getTrueValue()*
    std::pow(10.0,m_spline->f(log10X))*
    dlog10FluxdP*std::log(10)*m_parameter[iparam].getScale();
}


// double SplineFunction::
// integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;


optimizers::Function * SplineFunction::clone() const 
{
  return new SplineFunction(*this);
}

SplineFunction::SplineFunction(const SplineFunction& o):
  Function(o),
  m_initialized(o.m_initialized), m_npoints(o.m_npoints), m_cpp(o.m_cpp),
  m_cValues(o.m_cValues), m_spline(o.m_spline)
{
  if(m_spline)m_spline = new RegularSpline(*m_spline);
}

void SplineFunction::updateCache() const
{
  if(!m_initialized)
    throw optimizers::Exception("SplineFunction::updateCache: "
				"Not initialized!");

  bool must_update=false;
  if((m_spline==0)||(m_cValues.size() != m_npoints+2))
     must_update=true;
  else
    for(unsigned ipar=0;ipar<m_npoints+2;ipar++)
      if(m_parameter[ipar+unsigned(Index1)].getValue() != m_cValues[ipar])
	{
	  must_update=true;
	  break;
	}
  
  if(must_update)
    {
      m_cValues.resize(m_npoints+2);
      delete m_spline;
      std::vector<double> y(m_npoints);
      for(unsigned ipar=0;ipar<m_npoints+2;ipar++)
	m_cValues[ipar] = m_parameter[ipar+unsigned(Index1)].getValue();
      for(unsigned ipar=0;ipar<m_npoints;ipar++)
	y[ipar] = m_parameter[ipar+unsigned(Log10Flux0)].getTrueValue();
      double xmin = log10(m_parameter[LowerLimit].getTrueValue());
      double xmax = log10(m_parameter[UpperLimit].getTrueValue());
      if(m_cpp == CPP_CENTER)
	{
	  double dx2 = 0.5*(xmax-xmin)/double(m_npoints);
	  xmin += dx2;
	  xmax -= dx2;
	}
      double gamma1 = m_parameter[Index1].getTrueValue();
      double gammaN = m_parameter[IndexN].getTrueValue();
      m_spline = new RegularSpline(xmin, xmax, y,
				   RegularSpline::BCS_CLAMPED, gamma1,
				   RegularSpline::BCS_CLAMPED, gammaN,
				   true);
    }
}

void SplineFunction::setAttributes(const Attributes& attributes) 
{
  Function::setAttributes(attributes);

  Attributes::const_iterator iatt;

  iatt = attributes.find("npoints");
  if(iatt == attributes.end())
    {
      std::ostringstream str;
      str << "SplineFunction::setAttributes: attribute \"npoints\" not found.";
      throw optimizers::Exception(str.str());
    }

  m_npoints = ::strtoul(iatt->second.c_str(),0,10);
  if(m_npoints<2)
    {
      std::ostringstream str;
      str << "SplineFunction::setAttributes: must have >=2 parameters. "
	  << "Requested: " << m_npoints;
      throw optimizers::Exception(str.str());
    }

  iatt = attributes.find("control_points");
  if(iatt != attributes.end())
    {
      if(iatt->second == "edge")m_cpp=CPP_EDGE;
      else if(iatt->second == "center")m_cpp=CPP_CENTER;
      else
	{
	  std::ostringstream str;
	  str << "SplineFunction::setAttributes: "
	      << "unknown control point scheme: \""
	      << iatt->second << "\" (must be \"center\" or \"edge\")";
	  throw optimizers::Exception(str.str());
	}
    }

  typedef optimizers::Parameter Parameter;

  if(m_initialized)
    throw optimizers::Exception("SplineFunction::setAttributes: "
				"Already initialized!");

  m_initialized = true;

  setMaxNumParams(3+m_npoints+2);

  // unsigned decades=0;
  // for(unsigned i=m_npoints-1;i;i/=10)decades++;
  
  addParam(Parameter("Index1", -2.0, -50, 50, true, 0));
  addParam(Parameter("IndexN", std::numeric_limits<double>::infinity(), 
		     -50, 50, true));
  for(unsigned iparam=0;iparam<m_npoints;iparam++)
    {
      char buffer[80];
      sprintf(buffer, "Log10Flx%d", iparam);
      addParam(Parameter(buffer, -9, -100, 100, true));
    }
}

void SplineFunction::initializeParameters()
{
  Function::initializeParameters();

  if(!m_initialized)
    throw optimizers::Exception("SplineFunction::postInitialize: "
				"Not initialized!");
  
  double indexN = m_parameter[IndexN].getTrueValue();
  std::pair<double, double> indexNBds = m_parameter[IndexN].getBounds();
  if((indexN<indexNBds.first)||(indexN>indexNBds.second))
    {
      double gamma = m_parameter[Index1].getTrueValue();
      m_parameter[IndexN].setTrueValue(gamma);
      double dx = (std::log10(m_parameter[UpperLimit].getTrueValue()) -
		   std::log10(m_parameter[LowerLimit].getTrueValue()));
      if(m_cpp == CPP_EDGE)
	dx /= double(m_npoints-1);
      else
	dx /= double(m_npoints);

      double y0 = m_parameter[Log10Flux0].getTrueValue();
      for(unsigned iparam=1;iparam<m_npoints;iparam++)
	{
	  unsigned idx = unsigned(Log10Flux0)+iparam;
	  m_parameter[idx].setTrueValue(y0+gamma*double(iparam)*dx);
	}
    }
}

void SplineFunction::getAttributes(Attributes& attributes) const
{
  Function::getAttributes(attributes);

  std::ostringstream str;
  str << m_npoints;
  attributes["npoints"] = str.str();

  if(m_cpp == CPP_EDGE)
    attributes["control_points"] = "edge";
  else
    attributes["control_points"] = "center";
}

}
