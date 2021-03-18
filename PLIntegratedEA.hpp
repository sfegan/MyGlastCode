//-*-mode:c++; mode:font-lock;-*-

#ifndef PLINTEGRATEDEA_HPP
#define PLINTEGRATEDEA_HPP

#include <cmath>
#include <vector>

#include <irfInterface/IEfficiencyFactor.h>
#include <irfInterface/IrfsFactory.h>
#include <irfLoader/Loader.h>

#include "FunctionXY.hpp"

class PLIntegratedEA
{
  // This class helps manage the multitude of effective area integral
  // that are needed to handle the livetime and phi angle corrections
public:
  PLIntegratedEA(irfInterface::Irfs* irf,
		 double emin, double emax, double gamma,
		 double ebins_per_decade = 100.0);
  PLIntegratedEA(const std::vector<irfInterface::Irfs *>& irfs,
		 double emin, double emax, double gamma,
		 double ebins_per_decade = 100.0);
  ~PLIntegratedEA();
  void setROIConstantRadius(double radius);
  void setROIRadiusFunction(FunctionXY& radius_fn);
  void setROIPerIRFRadiusFunction(std::vector<FunctionXY*>& radius_fns);
  void setROIConstantPSFFraction(double fraction);
  double value(double costheta, double phi, double livetime_frac) const;
private:
  const std::vector<irfInterface::Irfs *>  m_irfs;
  double                                   m_lnemin;
  double                                   m_lnemax;
  double                                   m_dlne;
  double                                   m_one_plus_gamma;
  std::vector<FunctionXY*>                 m_radius_vs_log10e;
  double                                   m_psf_factor;
  double                                   m_pl_normalization;

  double swf(double lne) const { return std::exp(lne*m_one_plus_gamma); }
};

#endif // #ifndef PLINTEGRATEDEA_HPP
