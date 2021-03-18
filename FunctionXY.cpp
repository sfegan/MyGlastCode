//-*-mode:c++; mode:font-lock;-*-

#include <algorithm>
#include <stdexcept>

#include "FunctionXY.hpp"

FunctionXY::~FunctionXY()
{
  // nothing to see here
}

FunctionXYConst::~FunctionXYConst()
{
  // nothing to see here
}

double FunctionXYConst::f(double x)
{
  return m_y;
}

FunctionXY* FunctionXYConst::copy() const
{
  return new FunctionXYConst(*this);
}

FunctionXYInterpolate::
FunctionXYInterpolate(const std::vector<XYpoint>& xy, bool extrapolate):
  FunctionXY(), m_xy(xy), m_extrapolate(extrapolate) 
{ 
  std::sort(m_xy.begin(),m_xy.end()); 
}

FunctionXYInterpolate::~FunctionXYInterpolate()
{
  // nothing to see here
}

double FunctionXYInterpolate::f(double x)
{
  std::vector<XYpoint>::iterator i = 
    std::upper_bound(m_xy.begin(), m_xy.end(), XYpoint(x,0));
  if(i == m_xy.begin())
    {
      if(m_extrapolate)return m_xy.front().second;
      else throw std::range_error("FunctionXYInterpolate::f - value before beginning of range");
    }
  if(i == m_xy.end())
    {
      if(m_extrapolate)return m_xy.back().second;
      else throw std::range_error("FunctionXYInterpolate::f - value beyond end of range");
    }
  
  double xf = (i->first-x)/(i->first-(i-1)->first);
  return (i-1)->second*xf + i->second*(1-xf);
}

FunctionXY* FunctionXYInterpolate::copy() const
{
  return new FunctionXYInterpolate(*this);
}
