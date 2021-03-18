//-*-mode:c++; mode:font-lock;-*-

#include "Accumulator.hpp"
#include "PLIntegratedEA.hpp"

PLIntegratedEA::
PLIntegratedEA(const IRFs* irfs, double dcostheta, double costheta0, 
	       double log10emin, double log10emax, double gamma,
	       double roi_radius, double dlog10e):
  m_irfs(irfs), m_dct(dcostheta), m_ct0(costheta0), m_ncostheta(), m_cpt()
{
  std::vector<double> costheta;
  for(double ctheta=costheta0;ctheta<1+dcostheta*0.5;
      ctheta+=dcostheta)costheta.push_back(ctheta);
  m_ncostheta = costheta.size();
  m_cpt.resize(ncpt(),std::vector<double>(m_ncostheta));
  
  Accumulator a_norm;
  for(double loge=log10emin;loge<log10emax;loge+=dlog10e)
    {
      const double spec = std::pow(10.0,-(gamma-1)*loge);
      a_norm.add(spec);
    }
  
  for(unsigned icpt=0;icpt<ncpt();icpt++)
    {
      std::vector<double> cpt(m_ncostheta);
      for(unsigned ict=0;ict<m_ncostheta;ict++)
	{
	  const double ct = costheta[ict];
	  Accumulator a_int;
	  for(double loge=log10emin;loge<log10emax;loge+=dlog10e)
	    {
	      const double spec = std::pow(10.0,-(gamma-1)*loge);
	      const double wgt  = weight(loge,icpt);
	      const double area = m_irfs->ea()->value(loge,ct);
	      const double psf  = (roi_radius<M_PI)?
		m_irfs->psf()->integral(roi_radius,loge,ct):1.0;
	      a_int.add(spec*area*wgt*psf);
	    }
	  cpt[ict]=a_int.sum()/a_norm.sum();
	}
      m_cpt[icpt] = cpt;
    }
}

double PLIntegratedEA::weight(double log10e, unsigned icpt) const
{
  double w = 1.0;
  if(m_irfs->eff())
    {
      if(icpt%2 == 0)w = m_irfs->eff()->p1(log10e);
      else w = m_irfs->eff()->p0(log10e);
    }
  return w;
}

double PLIntegratedEA::
value(double costheta, double phi, double livetime_frac) const
{
  if(costheta<m_ct0)return 0;
  unsigned ict=lrint((costheta-m_ct0)/m_dct-0.5);
  double dct = costheta-m_ct0-double(ict)*m_dct;

  unsigned nc = ncpt();
  std::vector<double> x(nc,1.0);
  if(m_irfs->eff())
    {
      x[1] *= livetime_frac;
      if(m_irfs->phi())x[3] *= livetime_frac;
    }
  if(m_irfs->phi())
    {
      assert(0);
      double cp = cos(phi);
      double sp = sin(phi);
      x[0] *= cp;
      if(m_irfs->eff())x[1] *= cp, x[2] *= sp, x[3] *= sp;
      else x[1] *= sp;
    }

  Accumulator aeff;
  for(unsigned ic=0;ic<nc;ic++)
    {
      double f;
      if(ict>=(m_ncostheta-1))
	f = m_cpt[ic][m_ncostheta-1];
      else
	f = (m_cpt[ic][ict]*(1.0-dct)+m_cpt[ic][ict+1]*dct);
      aeff.add(f*x[ic]);
    }

  return aeff.sum();
}

