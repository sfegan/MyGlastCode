//-*-mode:c++; mode:font-lock;-*-

#ifndef PLINTEGRATEDEA_HPP
#define PLINTEGRATEDEA_HPP

#include <vector>

#include "IRF.hpp"

class PLIntegratedEA
{
  // This class helps manage the multitude of effective area integral
  // that are needed to handle the livetime and phi angle corrections
public:
  PLIntegratedEA(const IRFs* irfs, double dcostheta, double costheta0, 
		 double log10emin, double log10emax, double gamma,
		 double roi_radius = M_PI, double dlog10e = 0.01);

  unsigned ncpt() const 
  { 
    return (m_irfs->eff()?2:1)*(m_irfs->phi()?2:1); 
  }

  double weight(double log10e, unsigned icpt) const;
  double value(double costheta, double phi, double livetime_frac) const;
private:
  const IRFs* m_irfs;
  double m_dct;
  double m_ct0; 
  unsigned m_ncostheta;
  std::vector<std::vector<double> > m_cpt;
};

#endif // #ifndef PLINTEGRATEDEA_HPP
