//-*-mode:c++; mode:font-lock;-*-

/*   
   FunctionXY.hpp - Stephen Fegan 
                  - sfegan@llr.in2p3.fr
                  - 18 May 2011

   Class to encapsulate different 1-dimensional functions

   $Id: RegularSpline.hpp 1952 2010-07-08 07:17:12Z sfegan $
*/

#ifndef FUNCTION_XY_HPP
#define FUNCTION_XY_HPP

#include <vector>

class FunctionXY
{
public:	
  virtual ~FunctionXY();
  virtual double f(double x) = 0;
  virtual FunctionXY* copy() const = 0;
};

class FunctionXYConst: public FunctionXY
{
public:	
  FunctionXYConst(double y): FunctionXY(), m_y(y) { }
  virtual ~FunctionXYConst();
  virtual double f(double x);
  virtual FunctionXY* copy() const;
private:
  double m_y;
};

class FunctionXYInterpolate: public FunctionXY
{
public:	
  typedef std::pair<double,double> XYpoint;
  FunctionXYInterpolate(const std::vector<XYpoint>& xy, bool extrapolate=false);
  virtual ~FunctionXYInterpolate();
  virtual double f(double x);
  virtual FunctionXY* copy() const;
private:
  std::vector<XYpoint> m_xy;
  bool                 m_extrapolate;
};

#endif
