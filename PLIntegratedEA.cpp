//-*-mode:c++; mode:font-lock;-*-

#include <cmath>
#include <cassert>

#include "Accumulator.hpp"
#include "PLIntegratedEA.hpp"

PLIntegratedEA::
PLIntegratedEA(irfInterface::Irfs* irf,
	       double emin, double emax, double gamma,
	       double ebins_per_decade):
  m_irfs(1,irf), m_lnemin(std::log(emin)), m_lnemax(std::log(emax)), 
  m_dlne(std::log(10)/ebins_per_decade), m_one_plus_gamma(1.0-gamma),
  m_radius_vs_log10e(), m_psf_factor(1.0), m_pl_normalization()
{
  Accumulator a;
  for(double lne=m_lnemin+0.5*m_dlne;lne<m_lnemax;lne+=m_dlne)a.add(swf(lne));
  m_pl_normalization = a.sum();
}

PLIntegratedEA::
PLIntegratedEA(const std::vector<irfInterface::Irfs *>& irfs,
	       double emin, double emax, double gamma,
	       double ebins_per_decade):
  m_irfs(irfs), m_lnemin(std::log(emin)), m_lnemax(std::log(emax)), 
  m_dlne(std::log(10)/ebins_per_decade), m_one_plus_gamma(1.0-gamma),
  m_radius_vs_log10e(), m_psf_factor(1.0), m_pl_normalization()
{
  Accumulator a;
  for(double lne=m_lnemin+0.5*m_dlne;lne<m_lnemax;lne+=m_dlne)a.add(swf(lne));
  m_pl_normalization = a.sum();
}

PLIntegratedEA::~PLIntegratedEA()
{
  for(unsigned ifn=0;ifn<m_radius_vs_log10e.size();ifn++)
    delete m_radius_vs_log10e[ifn];
}

void PLIntegratedEA::setROIConstantRadius(double radius)
{
  for(unsigned ifn=0;ifn<m_radius_vs_log10e.size();ifn++)
    delete m_radius_vs_log10e[ifn];
  m_radius_vs_log10e.resize(m_irfs.size());
  for(unsigned iirf=0;iirf!=m_irfs.size(); iirf++)
    m_radius_vs_log10e[iirf] = new FunctionXYConst(radius);
  m_psf_factor = 1;
}

void PLIntegratedEA::setROIRadiusFunction(FunctionXY& radius_fn)
{
  for(unsigned ifn=0;ifn<m_radius_vs_log10e.size();ifn++)
    delete m_radius_vs_log10e[ifn];
  m_radius_vs_log10e.resize(m_irfs.size());
  for(unsigned iirf=0;iirf!=m_irfs.size(); iirf++)
    m_radius_vs_log10e[iirf] = radius_fn.copy();
  m_psf_factor = 1;
}

void PLIntegratedEA::
setROIPerIRFRadiusFunction(std::vector<FunctionXY*>& radius_fns)
{
  assert(radius_fns.size() == m_irfs.size());
  for(unsigned ifn=0;ifn<m_radius_vs_log10e.size();ifn++)
    delete m_radius_vs_log10e[ifn];
  m_radius_vs_log10e.resize(m_irfs.size());
  for(unsigned iirf=0;iirf!=m_irfs.size(); iirf++)
    m_radius_vs_log10e[iirf] = radius_fns[iirf]->copy();
  m_psf_factor = 1;
}

void PLIntegratedEA::setROIConstantPSFFraction(double fraction)
{
  for(unsigned ifn=0;ifn<m_radius_vs_log10e.size();ifn++)
    delete m_radius_vs_log10e[ifn];
  m_radius_vs_log10e.clear();
  m_psf_factor = fraction;  
}

double PLIntegratedEA::
value(double costheta, double phi, double livetime_frac) const
{
  const double phi_deg = phi*180.0/M_PI;
  const double theta_deg = std::acos(costheta)*180.0/M_PI;
  Accumulator swea;
  for(double lne=m_lnemin+0.5*m_dlne;lne<m_lnemax;lne+=m_dlne)
    {
      double energy = exp(lne);
      double ea(0);
      for (size_t iirf = 0; iirf < m_irfs.size(); iirf++) 
	{
	  double aperture(m_psf_factor);
	  if (!m_radius_vs_log10e.empty())
	    {
	      double radius = m_radius_vs_log10e[iirf]->f(lne/std::log(10));
	      if (radius < 180.0) 
		{
		  irfInterface::IPsf* psf = m_irfs[iirf]->psf();
		  aperture = 
		    psf->angularIntegral(energy, theta_deg, phi_deg, radius);
		}
	    }

	  const irfInterface::IEfficiencyFactor* eff = 
	    m_irfs[iirf]->efficiencyFactor();
	  double efficiency(1);
	  if(eff)efficiency = eff->value(energy, livetime_frac);

	  irfInterface::IAeff * aeff = m_irfs[iirf]->aeff();
	  ea += 
	    aeff->value(energy, theta_deg, phi_deg, 0 /* time unused */)
	    *aperture
	    *efficiency
	    *0.0001; // ST returns EA in cm^2
	}
      swea.add(swf(lne) * ea);
    }
  return swea.sum()/m_pl_normalization;
}
