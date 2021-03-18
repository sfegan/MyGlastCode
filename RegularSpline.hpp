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

   $Id: RegularSpline.hpp 1952 2010-07-08 07:17:12Z sfegan $

*/

#ifndef REGULAR_SPLINE_HPP
#define REGULAR_SPLINE_HPP

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

#endif // def REGULAR_SPLINE_HPP
