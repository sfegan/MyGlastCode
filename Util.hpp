//-*-mode:c++; mode:font-lock;-*-

#ifndef UTIL_HPP
#define UTIL_HPP

#include <utility>
#include <cmath>

inline double r2d(const double x)
{
  return x/M_PI*180.0;
}

inline double d2r(const double x)
{
  return x*M_PI/180.0;
}

inline double sqr(const double x)
{
  return x*x;
}

// Fast function which uses no trignometric functions to return a
// "spherical distance"-like value that can be used for sorting
// distances (only!)

inline double approx_angle(const double theta)
{
  const double z = std::cos(theta);
  const double rho = std::sin(theta);
  if(z>=0)
    if(z>rho)return rho/z;
    else return 2.0-z/rho;
  else
    {
      const double zz = -z;
      if(zz>rho)return 2.0+rho/zz;
      else return 4.0-z/rho;
    }
}  

inline double approx_sphere_dist_nt(const double sinra1,  const double cosra1, 
				    const double sindec1, const double cosdec1,
				    const double sinra2,  const double cosra2, 
				    const double sindec2, const double cosdec2)
{
  const double sindphi=sinra1*cosra2-cosra1*sinra2;
  const double cosdphi=cosra1*cosra2+sinra1*sinra2;
  const double sinth1=cosdec1;
  const double costh1=sindec1;
  const double sinth2=cosdec2;
  const double costh2=sindec2;
  const double x=sinth2*costh1-costh2*sinth1*cosdphi;
  const double y=sinth1*sindphi;
  const double z=costh2*costh1+sinth2*sinth1*cosdphi;
  const double rho = std::sqrt(x*x+y*y);
  if(z>=0)
    if(z>rho)return rho/z;
    else return 2.0-z/rho;
  else
    {
      const double zz = -z;
      if(zz>rho)return 2.0+rho/zz;
      else return 4.0-z/rho;
    }
}

inline double sphere_dist_nt(const double sinra1,  const double cosra1, 
			     const double sindec1, const double cosdec1,
			     const double sinra2,  const double cosra2, 
			     const double sindec2, const double cosdec2)
{
  const double sindphi=sinra1*cosra2-cosra1*sinra2;
  const double cosdphi=cosra1*cosra2+sinra1*sinra2;
  const double sinth1=cosdec1;
  const double costh1=sindec1;
  const double sinth2=cosdec2;
  const double costh2=sindec2;
  const double x=sinth2*costh1-costh2*sinth1*cosdphi;
  const double y=sinth1*sindphi;
  const double z=costh2*costh1+sinth2*sinth1*cosdphi;
  return std::atan2(std::sqrt(x*x+y*y),z);
}

inline double sphere_dist(const double ra1, const double dec1, 
			  const double ra2, const double dec2)
{
  const double sindphi=std::sin(ra1-ra2);
  const double cosdphi=std::cos(ra1-ra2);
  const double sinth1=std::cos(dec1);
  const double costh1=std::sin(dec1);
  const double sinth2=std::cos(dec2);
  const double costh2=std::sin(dec2);
  const double x=sinth2*costh1-costh2*sinth1*cosdphi;
  const double y=sinth1*sindphi;
  const double z=costh2*costh1+sinth2*sinth1*cosdphi;
  return std::atan2(std::sqrt(x*x+y*y),z);
}

template<typename T> inline void setCut(std::pair<bool,T>& cut, const T& x)
{
  cut.first = true;
  cut.second = x;
}

#endif
