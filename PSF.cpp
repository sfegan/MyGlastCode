#include <cmath>
#include <iostream>
#include <sstream>

#include "FITS.hpp"
#include "Util.hpp"
#include "PSF.hpp"
#include "Accumulator.hpp"

PSFCalc::~PSFCalc()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// 
// ISOTROPIC BACKGROUND
// 
// ----------------------------------------------------------------------------

Like_IsotropicBkg::
Like_IsotropicBkg(const std::vector<double>& d, double dmax)
  : PSFCalc(), m_nd(d.size()), m_dmax2(dmax*dmax)
{
  // nothing to see here
}

Like_IsotropicBkg::~Like_IsotropicBkg()
{
  // nothing to see here
}

unsigned Like_IsotropicBkg::numDim()
{
  return 1;
}

double Like_IsotropicBkg::Pdiff(const std::vector<double>& x,
				const std::vector<double>& p)
{
  const double N  = double(m_nd)*p[0];
  const double norm_b = N/(M_PI*m_dmax2);
  return norm_b;
}

double Like_IsotropicBkg::Pint(const std::vector<double>& x,
			       const std::vector<double>& x0,
			       const std::vector<double>& p)
{
  const double N  = double(m_nd)*p[0];
  const double norm_b = N/m_dmax2;
  return norm_b*(x[0]*x[0] - x0[0]*x0[0]);
}

std::string Like_IsotropicBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_IsotropicBkg: unknown dimension");
    }
}

bool Like_IsotropicBkg::
canCalcDFDP()
{
  return false;
}

double Like_IsotropicBkg::
f(const std::vector<double>& p)
{
  const double N  = double(m_nd)*p[0];
  const double norm_b = N/(M_PI*m_dmax2);
  double log_like = double(m_nd)*log(norm_b) - N;
  return -log_like;
}

double Like_IsotropicBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_IsotropicBkg: dfdp called");
}

unsigned Like_IsotropicBkg::numParam()
{
  return 1;
}

std::string Like_IsotropicBkg::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "N";
    default: throw std::string("Like_IsotropicBkg: unknown parameter");
    }
}

double Like_IsotropicBkg::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 1.0;
    default: throw std::string("Like_IsotropicBkg: unknown parameter");
    }  
}

double Like_IsotropicBkg::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.5;
    default: throw std::string("Like_IsotropicBkg: unknown parameter");
    }  
}

double Like_IsotropicBkg::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 2.0;
    default: throw std::string("Like_IsotropicBkg: unknown parameter");
    }  
}

double Like_IsotropicBkg::
modelRP(const std::vector<double>& p, const double P)
{
  return 0;
}

void Like_IsotropicBkg::
modelRP(double& rp, double& rp_err,
       const std::vector<double>& p,
       const std::vector<double>& p_err,
       const std::vector<std::vector<double> > p_cov,
       const double P)
{
  rp = rp_err = 0;
  return;
}

// ----------------------------------------------------------------------------
// 
// GAUSSIAN
// 
// ----------------------------------------------------------------------------

Like_GaussianPSFWithBkg::
Like_GaussianPSFWithBkg(const std::vector<double>& d, double dmax)
  : PSFCalc(), m_d2(), m_dmax2(dmax*dmax)
{
  const unsigned nd = d.size();
  m_d2.resize(nd);
  for(unsigned id=0;id<nd;id++)m_d2[id] = d[id]*d[id];
}

Like_GaussianPSFWithBkg::~Like_GaussianPSFWithBkg()
{
  // nothing to see here
}

unsigned Like_GaussianPSFWithBkg::numDim()
{
  return 1;
}

double Like_GaussianPSFWithBkg::Pdiff(const std::vector<double>& x,
				      const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma = p[0];
  const double minus_one_over_2_sigma_sq = -0.5/(sigma*sigma);
  
  const double Ps = p[1];
  const double Pb = 1-Ps;
  const double N  = double(nd)*p[2];

  const double norm_s = 
    N*Ps/((2.0*M_PI*sigma*sigma)
	  *(1-std::exp(m_dmax2*minus_one_over_2_sigma_sq)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  return norm_s*std::exp(x[0]*x[0]*minus_one_over_2_sigma_sq) + norm_b;
}

double Like_GaussianPSFWithBkg::Pint(const std::vector<double>& x,
				     const std::vector<double>& x0,
				     const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma = p[0];
  const double minus_one_over_2_sigma_sq = -0.5/(sigma*sigma);
  
  const double Ps = p[1];
  const double Pb = 1-Ps;
  const double N  = double(nd)*p[2];

  const double norm_s = 
    N*Ps/(1-std::exp(m_dmax2*minus_one_over_2_sigma_sq));
  const double norm_b = 
    N*Pb/m_dmax2;

  return norm_s*(std::exp(x0[0]*x0[0]*minus_one_over_2_sigma_sq)
		 - std::exp(x[0]*x[0]*minus_one_over_2_sigma_sq))
    + norm_b*(x[0]*x[0] - x0[0]*x0[0]);
}

std::string Like_GaussianPSFWithBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_GaussianPSFWithBkg: unknown dimension");
    }
}

bool Like_GaussianPSFWithBkg::
canCalcDFDP()
{
  return false;
}

double Like_GaussianPSFWithBkg::
f(const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma = p[0];
  const double minus_one_over_2_sigma_sq = -0.5/(sigma*sigma);
  
  const double Ps = p[1];
  const double Pb = 1-Ps;
  const double N  = double(nd)*p[2];

  const double norm_s = 
    N*Ps/((2.0*M_PI*sigma*sigma)
	  *(1-std::exp(m_dmax2*minus_one_over_2_sigma_sq)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  double log_like = 0;

  for(unsigned id=0;id<nd;id++)
    log_like += log(norm_s*std::exp(m_d2[id]*minus_one_over_2_sigma_sq) 
		    + norm_b);
  log_like -= N;

  return -log_like;
}

double Like_GaussianPSFWithBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_GaussianPSFWithBkg: dfdp called");
}

unsigned Like_GaussianPSFWithBkg::numParam()
{
  return 3;
}

std::string Like_GaussianPSFWithBkg::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "sigma";
    case 1: return "P_s";
    case 2: return "N";
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }
}

double Like_GaussianPSFWithBkg::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.2;
    case 1: return 0.5;
    case 2: return 1.0;
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }  
}

double Like_GaussianPSFWithBkg::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.0;
    case 1: return 0.0;
    case 2: return 0.5;
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }  
}

double Like_GaussianPSFWithBkg::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 5.0;
    case 1: return 2.0;
    case 2: return 2.0;
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }  
}

double Like_GaussianPSFWithBkg::
modelRP(const std::vector<double>& p, const double P)
{
  return p[0]*std::sqrt(-2*log(1-P));
}

void Like_GaussianPSFWithBkg::
modelRP(double& rp, double& rp_err,
       const std::vector<double>& p,
       const std::vector<double>& p_err,
       const std::vector<std::vector<double> > p_cov,
       const double P)
{
  const double F = std::sqrt(-2*log(1-P));
  rp = p[0]*F;
  rp_err = p_err[0]*F;
}

// ----------------------------------------------------------------------------
// 
// DOUBLE GAUSSIAN
// 
// ----------------------------------------------------------------------------

Like_DblGaussianPSFWithBkg::
Like_DblGaussianPSFWithBkg(const std::vector<double>& d, double dmax)
  : PSFCalc(), m_d2(), m_dmax2(dmax*dmax)
{
  const unsigned nd = d.size();
  m_d2.resize(nd);
  for(unsigned id=0;id<nd;id++)m_d2[id] = d[id]*d[id];
}

Like_DblGaussianPSFWithBkg::~Like_DblGaussianPSFWithBkg()
{
  // nothing to see here
}

unsigned Like_DblGaussianPSFWithBkg::numDim()
{
  return 1;
}

std::string Like_DblGaussianPSFWithBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_GaussianPSFWithBkg: unknown dimension");
    }
}

bool Like_DblGaussianPSFWithBkg::
canCalcDFDP()
{
  return false;
}

double Like_DblGaussianPSFWithBkg::Pdiff(const std::vector<double>& x,
					 const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma1 = p[0];
  const double sigma2 = sigma1*p[1];
  const double minus_one_over_2_sigma1_sq = -0.5/(sigma1*sigma1);
  const double minus_one_over_2_sigma2_sq = -0.5/(sigma2*sigma2);

  const double Ps  = p[2];
  const double Pb  = 1-Ps;
  const double Ps1 = Ps*p[3];
  const double Ps2 = Ps-Ps1;
  const double N   = double(nd)*p[4];

  const double norm_s1 = 
    N*Ps1/((2*M_PI*sigma1*sigma1)
	   *(1-std::exp(m_dmax2*minus_one_over_2_sigma1_sq)));
  const double norm_s2 = 
    N*Ps2/((2*M_PI*sigma2*sigma2)
	   *(1-std::exp(m_dmax2*minus_one_over_2_sigma2_sq)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  return norm_s1*std::exp(x[0]*x[0]*minus_one_over_2_sigma1_sq) 
       + norm_s2*std::exp(x[0]*x[0]*minus_one_over_2_sigma2_sq) + norm_b;
}

double Like_DblGaussianPSFWithBkg::Pint(const std::vector<double>& x,
					const std::vector<double>& x0,
					const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma1 = p[0];
  const double sigma2 = sigma1*p[1];
  const double minus_one_over_2_sigma1_sq = -0.5/(sigma1*sigma1);
  const double minus_one_over_2_sigma2_sq = -0.5/(sigma2*sigma2);
  
  const double Ps  = p[2];
  const double Pb  = 1-Ps;
  const double Ps1 = Ps*p[3];
  const double Ps2 = Ps-Ps1;
  const double N   = double(nd)*p[4];

  const double norm_s1 = 
    N*Ps1/(1-std::exp(m_dmax2*minus_one_over_2_sigma1_sq));
  const double norm_s2 = 
    N*Ps2/(1-std::exp(m_dmax2*minus_one_over_2_sigma2_sq));
  const double norm_b = 
    N*Pb/m_dmax2;

  return norm_s1*(std::exp(x0[0]*x0[0]*minus_one_over_2_sigma1_sq)
		  - std::exp(x[0]*x[0]*minus_one_over_2_sigma1_sq))
    + norm_s2*(std::exp(x0[0]*x0[0]*minus_one_over_2_sigma2_sq)
	       - std::exp(x[0]*x[0]*minus_one_over_2_sigma2_sq))
    + norm_b*(x[0]*x[0] - x0[0]*x0[0]);
}

double Like_DblGaussianPSFWithBkg::
f(const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma1 = p[0];
  const double sigma2 = sigma1*p[1];
  const double minus_one_over_2_sigma1_sq = -0.5/(sigma1*sigma1);
  const double minus_one_over_2_sigma2_sq = -0.5/(sigma2*sigma2);
  
  const double Ps  = p[2];
  const double Pb  = 1-Ps;
  const double Ps1 = Ps*p[3];
  const double Ps2 = Ps-Ps1;
  const double N   = double(nd)*p[4];

  const double norm_s1 = 
    N*Ps1/((2*M_PI*sigma1*sigma1)
	   *(1-std::exp(m_dmax2*minus_one_over_2_sigma1_sq)));
  const double norm_s2 = 
    N*Ps2/((2*M_PI*sigma2*sigma2)
	   *(1-std::exp(m_dmax2*minus_one_over_2_sigma2_sq)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  double log_like = 0;

  for(unsigned id=0;id<nd;id++)
    log_like += log(norm_s1*std::exp(m_d2[id]*minus_one_over_2_sigma1_sq)
		    + norm_s2*std::exp(m_d2[id]*minus_one_over_2_sigma2_sq)
		    + norm_b);
  log_like -= N;

  return -log_like;

}

double Like_DblGaussianPSFWithBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_GaussianPSFWithBkg: dfdp called");
}

unsigned Like_DblGaussianPSFWithBkg::numParam()
{
  return 5;
}

std::string Like_DblGaussianPSFWithBkg::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "sig1";
    case 1: return "sig2/sig1";
    case 2: return "P_s";
    case 3: return "P_s2/Ps";
    case 4: return "N";
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }
}

double Like_DblGaussianPSFWithBkg::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.2;
    case 1: return 2.0;
    case 2: return 0.5;
    case 3: return 0.25;
    case 4: return 1.0;
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }  
}

double Like_DblGaussianPSFWithBkg::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.0;
    case 1: return 1.0;
    case 2: return 0.0;
    case 3: return 0.0;
    case 4: return 0.5;
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }  
}

double Like_DblGaussianPSFWithBkg::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 5.0;
    case 1: return 100.0;
    case 2: return 1.1;
    case 3: return 1.1;
    case 4: return 2.0;
    default: throw std::string("Like_GaussianPSFWithBkg: unknown parameter");
    }  
}

double Like_DblGaussianPSFWithBkg::
modelRP(const std::vector<double>& p, const double P)
{
  const double sigma1 = p[0];
  const double sigma2 = sigma1*p[1];
  const double minus_one_over_2_sigma1_sq = -0.5/(sigma1*sigma1);
  const double minus_one_over_2_sigma2_sq = -0.5/(sigma2*sigma2);
  
  const double Ps1 = p[3];
  const double Ps2 = 1-Ps1;

  double rP0 = 0;
  double rP1 = p[0]*std::sqrt(-2*log(1-P));
  double PrP = 1
    - Ps1*std::exp(rP1*rP1*minus_one_over_2_sigma1_sq)
    - Ps2*std::exp(rP1*rP1*minus_one_over_2_sigma2_sq);
  while(PrP < P)
    {
      rP0=rP1;
      rP1*=2.0;
      PrP = 1
	- Ps1*std::exp(rP1*rP1*minus_one_over_2_sigma1_sq)
	- Ps2*std::exp(rP1*rP1*minus_one_over_2_sigma2_sq);
    }
  while(fabs(rP0-rP1)>1e-8)
    {
      const double rPt = 0.5*(rP0+rP1);
      PrP = 1
	- Ps1*std::exp(rPt*rPt*minus_one_over_2_sigma1_sq)
	- Ps2*std::exp(rPt*rPt*minus_one_over_2_sigma2_sq);
      if(PrP > P)rP1=rPt;
      else rP0=rPt;
    }

  return 0.5*(rP0+rP1);
}

void Like_DblGaussianPSFWithBkg::
modelRP(double& rp, double& rp_err,
	const std::vector<double>& p,
	const std::vector<double>& p_err,
	const std::vector<std::vector<double> > p_cov,
	const double P)
{
  rp = modelRP(p,P);

  std::vector<double> pt(p);
  double dx = 0;

  dx = p[0]*1e-5;
  pt[0] += dx;
  double dr_dp0 = (modelRP(pt,P)-rp)/dx;
  pt[0] = p[0];

  dx = p[1]*1e-5;
  pt[1] += dx;
  double dr_dp1 = (modelRP(pt,P)-rp)/dx;
  pt[1] = p[1];

  dx = p[2]*1e-5;
  pt[2] += dx;
  double dr_dp2 = (modelRP(pt,P)-rp)/dx;
  pt[2] = p[2];

  rp_err = std::sqrt(p_cov[0][0]*dr_dp0*dr_dp0
		     + p_cov[1][1]*dr_dp1*dr_dp1
		     + p_cov[2][2]*dr_dp2*dr_dp2
		     + 2.0*p_cov[0][1]*dr_dp0*dr_dp1
		     + 2.0*p_cov[1][2]*dr_dp1*dr_dp2
		     + 2.0*p_cov[2][0]*dr_dp2*dr_dp0);
}

// ----------------------------------------------------------------------------
// 
// MULTI GAUSSIAN
// 
// ----------------------------------------------------------------------------

Like_MultiGaussianPSFWithBkg::
Like_MultiGaussianPSFWithBkg(unsigned n,
			     const std::vector<double>& d, double dmax)
  : PSFCalc(), m_n(n), m_d2(), m_dmax2(dmax*dmax)
{
  const unsigned nd = d.size();
  m_d2.resize(nd);
  for(unsigned id=0;id<nd;id++)m_d2[id] = d[id]*d[id];
}

Like_MultiGaussianPSFWithBkg::~Like_MultiGaussianPSFWithBkg()
{
  // nothing to see here
}

unsigned Like_MultiGaussianPSFWithBkg::numDim()
{
  return 1;
}

std::string Like_MultiGaussianPSFWithBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_GaussianPSFWithBkg: unknown dimension");
    }
}

bool Like_MultiGaussianPSFWithBkg::
canCalcDFDP()
{
  return false;
}

double Like_MultiGaussianPSFWithBkg::Pdiff(const std::vector<double>& x,
					 const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();
  const double x2 = x[0]*x[0];

  double sigmaC = 1.0;
  double PsC    = 0.0;
  if(m_n)PsC = p[1];

  double P = (1-PsC)/(M_PI*m_dmax2);
  for(unsigned i=0;i<m_n;i++)
    {
      const double sigmai  = sigmaC*p[i*2];
      sigmaC = sigmai;
      const double sigmai2 = sigmai*sigmai;
      double Psi           = PsC;
      if(i != m_n-1)Psi *= p[i*2+3], PsC -= Psi;
      const double m1_2si2 = -0.5/sigmai2;
      const double norm    = Psi/(2*M_PI*sigmai2*(1-std::exp(m_dmax2*m1_2si2)));
      P += norm*std::exp(x2*m1_2si2);
    }

  return P*double(nd)*p[2*m_n];
}

double Like_MultiGaussianPSFWithBkg::Pint(const std::vector<double>& x,
					const std::vector<double>& x0,
					const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();
  const double x2 = x[0]*x[0];
  const double x02 = x0[0]*x0[0];

  double sigmaC = 1.0;
  double PsC    = 0.0;
  if(m_n)PsC = p[1];

  double I = (1-PsC)*(x2 - x02)/m_dmax2;
  for(unsigned i=0;i<m_n;i++)
    {
      const double sigmai  = sigmaC*p[i*2];
      sigmaC = sigmai;
      const double sigmai2 = sigmai*sigmai;
      double Psi           = PsC;
      if(i != m_n-1)Psi *= p[i*2+3], PsC -= Psi;
      const double m1_2si2 = -0.5/sigmai2;
      const double norm    = Psi/(1-std::exp(m_dmax2*m1_2si2));
      I += norm*(std::exp(x02*m1_2si2) - std::exp(x2*m1_2si2));
    }

  return I*double(nd)*p[2*m_n];
}

double Like_MultiGaussianPSFWithBkg::
f(const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();
  const double N = double(nd)*p[2*m_n];

  double log_like = 0;
  for(unsigned id=0;id<nd;id++)
    {
      double sigmaC = 1.0;
      double PsC    = 0.0;
      if(m_n)PsC = p[1];

      double Li = (1-PsC)/(M_PI*m_dmax2);

      for(unsigned i=0;i<m_n;i++)
	{
	  const double sigmai  = sigmaC*p[i*2];
	  sigmaC = sigmai;
	  const double sigmai2 = sigmai*sigmai;
	  double Psi           = PsC;
	  if(i != m_n-1)Psi *= p[i*2+3], PsC -= Psi;
	  const double m1_2si2 = -0.5/sigmai2;
	  const double norm    = Psi/(2*M_PI*sigmai2*(1-std::exp(m_dmax2*m1_2si2)));
	  Li += norm*std::exp(m_d2[id]*m1_2si2);
	}

      log_like += log(Li);
    }

  log_like += double(nd)*log(N) - N;
  return -log_like;
}

double Like_MultiGaussianPSFWithBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_MultiGaussianPSFWithBkg: dfdp called");
}

unsigned Like_MultiGaussianPSFWithBkg::numParam()
{
  return 2*m_n+1;
}

std::string Like_MultiGaussianPSFWithBkg::pName(unsigned iparam)
{
  if(iparam<m_n*2)
    {
      if(iparam == 0)return "sig1";
      else if(iparam == 1)return "P_s";
      else if(iparam%2 == 0)
	{
	  unsigned ig = iparam/2 + 1;
	  std::ostringstream os;
	  os << "sig" << ig << "/sig" << ig-1;
	  return os.str();
	}
      else
	{
	  unsigned ig = iparam/2 + 1;
	  std::ostringstream os;
	  os << "P_s" << ig-1 << "/P_s" << ig;
	  return os.str();
	}
    }
  else if(iparam==m_n*2)return "N";
  else throw std::string("Like_MultiGaussianPSFWithBkg: unknown parameter");
}

double Like_MultiGaussianPSFWithBkg::p0(unsigned iparam)
{
  if(iparam<m_n*2)
    {
      if(iparam == 0)return 0.2;
      else if(iparam == 1)return 0.5;
      else if(iparam%2 == 0)return 2.0;
      else return 0.25;
    }
  else if(iparam==m_n*2)return 1.0;
  else throw std::string("Like_MultiGaussianPSFWithBkg: unknown parameter");
}  

double Like_MultiGaussianPSFWithBkg::plo(unsigned iparam)
{
  if(iparam<m_n*2)
    {
      if(iparam == 0)return 0.0;
      else if(iparam%2 == 0)return 1.0;
      else return 0.0;
    }
  else if(iparam==m_n*2)return 0.5;
  else throw std::string("Like_MultiGaussianPSFWithBkg: unknown parameter");
}

double Like_MultiGaussianPSFWithBkg::phi(unsigned iparam)
{
  if(iparam<m_n*2)
    {
      if(iparam == 0)return 5.0;
      else if(iparam%2 == 0)return 100.0;
      else return 1.1;
    }
  else if(iparam==m_n*2)return 1.1;
  else throw std::string("Like_MultiGaussianPSFWithBkg: unknown parameter");
}

double Like_MultiGaussianPSFWithBkg::
I(double x, const std::vector<double>& p) const
{
  const double x2 = x*x;

  double sigmaC = 1.0;
  double PsC    = 1.0;

  double I = 1.0;
  for(unsigned i=0;i<m_n;i++)
    {
      const double sigmai  = sigmaC*p[i*2];
      sigmaC = sigmai;
      const double sigmai2 = sigmai*sigmai;
      double Psi           = PsC;
      if(i != m_n-1)Psi *= p[i*2+3], PsC -= Psi;
      const double m1_2si2 = -0.5/sigmai2;
      I -= Psi*std::exp(x2*m1_2si2);
    }

  return I;
}


double Like_MultiGaussianPSFWithBkg::
modelRP(const std::vector<double>& p, const double P)
{
  double rP0 = 0;
  double rP1 = p[0]*std::sqrt(-2*log(1-P));
  double PrP = I(rP1,p);
  while(PrP < P)
    {
      rP0=rP1;
      rP1*=2.0;
      PrP = I(rP1,p);
    }
  while(fabs(rP0-rP1)>1e-8)
    {
      const double rPt = 0.5*(rP0+rP1);
      PrP = I(rPt,p);
      if(PrP > P)rP1=rPt;
      else rP0=rPt;
    }

  return 0.5*(rP0+rP1);
}

void Like_MultiGaussianPSFWithBkg::
modelRP(double& rp, double& rp_err,
	const std::vector<double>& p,
	const std::vector<double>& p_err,
	const std::vector<std::vector<double> > p_cov,
	const double P)
{
  rp = modelRP(p,P);

  std::vector<double> dr_dpi(m_n*2);

  for(unsigned i=0;i<m_n*2;i++)
    {
      std::vector<double> pt(p);
      double dx = 0;
      dx = p[i]*1e-5;
      pt[i] += dx;
      dr_dpi[i] = (modelRP(pt,P)-rp)/dx;
    }

  double rp2 = 0;
  for(unsigned i=0;i<m_n*2;i++)
    for(unsigned j=0;j<m_n*2;j++)
      rp2 += p_cov[i][j]*dr_dpi[i]*dr_dpi[j];
  rp_err = std::sqrt(rp2);
}

// ----------------------------------------------------------------------------
// 
// BURDETT FUNCTION
// 
// ----------------------------------------------------------------------------

Like_BurdettPSFWithBkg::
Like_BurdettPSFWithBkg(const std::vector<double>& d, double dmax)
  : PSFCalc(), m_d2(), m_dmax2(dmax*dmax)
{
  const unsigned nd = d.size();
  m_d2.resize(nd);
  for(unsigned id=0;id<nd;id++)m_d2[id] = d[id]*d[id];
}

Like_BurdettPSFWithBkg::~Like_BurdettPSFWithBkg()
{
  // nothing to see here
}

unsigned Like_BurdettPSFWithBkg::numDim()
{
  return 1;
}

double Like_BurdettPSFWithBkg::Pdiff(const std::vector<double>& x,
					const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma = p[0];
  const double gamma = p[1];

  const double Ps = p[2];
  const double Pb = 1.0-Ps;
  const double N  = double(nd)*p[3];

  const double one_over_2_gamma_sigma_sq = 0.5/(gamma*sigma*sigma);
  const double minus_gamma = -gamma;
  const double minus_gamma_minus_one = -(gamma-1.0);
  const double norm_s = 
    N*Ps*(1.0-1.0/gamma)
    / ((2.0*M_PI*sigma*sigma)
       *(1.0-std::pow(1.0+m_dmax2*one_over_2_gamma_sigma_sq,
				       minus_gamma_minus_one)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  return norm_s*pow(1.0+x[0]*x[0]*one_over_2_gamma_sigma_sq,
			       minus_gamma) + norm_b;
}

double Like_BurdettPSFWithBkg::Pint(const std::vector<double>& x,
				       const std::vector<double>& x0,
				       const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma = p[0];
  const double gamma = p[1];

  const double Ps = p[2];
  const double Pb = 1.0-Ps;
  const double N  = double(nd)*p[3];

  const double one_over_2_gamma_sigma_sq = 0.5/(gamma*sigma*sigma);
  const double minus_gamma_minus_one = -(gamma-1.0);
  const double norm_s = 
    N*Ps/(1.0-pow(1.0+m_dmax2*one_over_2_gamma_sigma_sq,
		   minus_gamma_minus_one));
  const double norm_b = 
    N*Pb/m_dmax2;

  return norm_s*(pow(1.0+x0[0]*x0[0]*one_over_2_gamma_sigma_sq,
		     minus_gamma_minus_one) - 
		 pow(1.0+x[0]*x[0]*one_over_2_gamma_sigma_sq,
		     minus_gamma_minus_one)) +
    norm_b*(x[0]*x[0] - x0[0]*x0[0]);
}

std::string Like_BurdettPSFWithBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_BurdettPSFWithBkg: unknown dimension");
    }
}

bool Like_BurdettPSFWithBkg::
canCalcDFDP()
{
  return false;
}

double Like_BurdettPSFWithBkg::
f(const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double sigma = p[0];
  const double gamma = p[1];

  const double Ps = p[2];
  const double Pb = 1.0-Ps;
  const double N  = double(nd)*p[3];

  const double one_over_2_gamma_sigma_sq = 0.5/(gamma*sigma*sigma);
  const double minus_gamma = -gamma;
  const double minus_gamma_minus_one = -(gamma-1.0);
  const double norm_s = 
    N*Ps*(1.0-1.0/gamma)
    / ((2.0*M_PI*sigma*sigma)*(1.0-pow(1.0+m_dmax2*one_over_2_gamma_sigma_sq,
				       minus_gamma_minus_one)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  double log_like = 0;

  for(unsigned id=0;id<nd;id++)
    log_like += log(norm_s*pow(1.0+m_d2[id]*one_over_2_gamma_sigma_sq,
			       minus_gamma) + norm_b);
  log_like -= N;

  return -log_like;
}

double Like_BurdettPSFWithBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_BurdettPSFWithBkg: dfdp called");
}

unsigned Like_BurdettPSFWithBkg::numParam()
{
  return 4;
}

std::string Like_BurdettPSFWithBkg::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "sigma";
    case 1: return "gamma";
    case 2: return "P_s";
    case 3: return "N";
    default: throw std::string("Like_BurdettPSFWithBkg: unknown parameter");
    }
}

double Like_BurdettPSFWithBkg::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.1;
    case 1: return 2.0;
    case 2: return 0.5;
    case 3: return 1.0;
    default: throw std::string("Like_BurdettPSFWithBkg: unknown parameter");
    }  
}

double Like_BurdettPSFWithBkg::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.0;
    case 1: return 1.0;
    case 2: return 0.0;
    case 3: return 0.5;
    default: throw std::string("Like_BurdettPSFWithBkg: unknown parameter");
    }  
}

double Like_BurdettPSFWithBkg::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 5.0;
    case 1: return 5.0;
    case 2: return 2.0;
    case 3: return 2.0;
    default: throw std::string("Like_BurdettPSFWithBkg: unknown parameter");
    }  
}

double Like_BurdettPSFWithBkg::
modelRP(const std::vector<double>& p, const double P)
{
  const double sigma = p[0];
  const double gamma = p[1];
  return sigma*std::sqrt(2.0*gamma*(pow(1.0-P,-1.0/(gamma-1.0))-1.0));
}

void Like_BurdettPSFWithBkg::
modelRP(double& rp, double& rp_err,
	const std::vector<double>& p,
	const std::vector<double>& p_err,
	const std::vector<std::vector<double> > p_cov,
	const double P)
{
  const double sigma = p[0];
  const double gamma = p[1];
  
  double gmo = gamma-1.0;
  double omP = 1.0-P;

  double powOmP = pow(omP,-1.0/gmo);

  double dr_dsigma = 
    std::sqrt(2.0*gamma*(powOmP-1.0));
  double dr_dgamma = 
    sigma/dr_dsigma*(powOmP*(1+gamma*log(omP)/(gmo*gmo))-1.0);

  rp = sigma * dr_dsigma;

  rp_err = std::sqrt(p_cov[0][0]*dr_dsigma*dr_dsigma
		     + p_cov[1][1]*dr_dgamma*dr_dgamma
		     + 2.0*p_cov[0][1]*dr_dsigma*dr_dgamma);
}

// ----------------------------------------------------------------------------
// 
// KING PROFILE
// 
// ----------------------------------------------------------------------------

Like_KingPSFWithBkg::
Like_KingPSFWithBkg(const std::vector<double>& d, double dmax,
		     double P)
  : PSFCalc(), m_P(P), m_d2(), m_dmax2(dmax*dmax)
{
  const unsigned nd = d.size();
  m_d2.resize(nd);
  for(unsigned id=0;id<nd;id++)m_d2[id] = d[id]*d[id];
}

Like_KingPSFWithBkg::~Like_KingPSFWithBkg()
{
  // nothing to see here
}

unsigned Like_KingPSFWithBkg::numDim()
{
  return 1;
}

double Like_KingPSFWithBkg::Pdiff(const std::vector<double>& x,
					const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double rPm   = p[0];
  const double gamma = 1.0/p[1];
  const double sigma = _sigma(rPm,gamma);

  const double Ps = p[2];
  const double Pb = 1.0-Ps;
  const double N  = double(nd)*p[3];

  const double sigma_sq = sigma*sigma;
  const double one_over_2_sigma_sq = 0.5/sigma_sq;
  const double minus_gamma = -gamma;
  const double minus_gamma_minus_one = -(gamma-1.0);
  const double norm_s = 
    N*Ps*(gamma-1.0)
    / ((2.0*M_PI*sigma_sq)*(1.0-pow(1.0+m_dmax2*one_over_2_sigma_sq,
				       minus_gamma_minus_one)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  return norm_s*pow(1.0+x[0]*x[0]*one_over_2_sigma_sq,
			       minus_gamma) + norm_b;
}

double Like_KingPSFWithBkg::Pint(const std::vector<double>& x,
				       const std::vector<double>& x0,
				       const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double rPm   = p[0];
  const double gamma = 1.0/p[1];
  const double sigma = _sigma(rPm,gamma);

  const double Ps = p[2];
  const double Pb = 1.0-Ps;
  const double N  = double(nd)*p[3];

  const double sigma_sq = sigma*sigma;
  const double one_over_2_sigma_sq = 0.5/sigma_sq;
  const double minus_gamma_minus_one = -(gamma-1.0);
  const double norm_s = 
    N*Ps/(1.0-pow(1.0+m_dmax2*one_over_2_sigma_sq,
		  minus_gamma_minus_one));
  const double norm_b = 
    N*Pb/m_dmax2;

  return norm_s*(pow(1.0+x0[0]*x0[0]*one_over_2_sigma_sq,
		     minus_gamma_minus_one) - 
		 pow(1.0+x[0]*x[0]*one_over_2_sigma_sq,
		     minus_gamma_minus_one)) +
    norm_b*(x[0]*x[0] - x0[0]*x0[0]);
}

std::string Like_KingPSFWithBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_KingPSFWithBkg: unknown dimension");
    }
}

bool Like_KingPSFWithBkg::
canCalcDFDP()
{
  return false;
}

double Like_KingPSFWithBkg::
f(const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double rPm   = p[0];
  const double gamma = 1.0/p[1];
  const double sigma = _sigma(rPm,gamma);

  const double Ps = p[2];
  const double Pb = 1.0-Ps;
  const double N  = double(nd)*p[3];

  const double sigma_sq = sigma*sigma;
  const double one_over_2_sigma_sq = 0.5/sigma_sq;
  const double minus_gamma = -gamma;
  const double minus_gamma_minus_one = -(gamma-1.0);
  const double norm_s = 
    N*Ps*(gamma-1.0)
    / ((2.0*M_PI*sigma_sq)*(1.0-pow(1.0+m_dmax2*one_over_2_sigma_sq,
				    minus_gamma_minus_one)));
  const double norm_b = 
    N*Pb/(M_PI*m_dmax2);

  double log_like = 0;

  for(unsigned id=0;id<nd;id++)
    log_like += log(norm_s*pow(1.0+m_d2[id]*one_over_2_sigma_sq,
			       minus_gamma) + norm_b);
  log_like -= N;

  return -log_like;
}

double Like_KingPSFWithBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_KingPSFWithBkg: dfdp called");
}

unsigned Like_KingPSFWithBkg::numParam()
{
  return 4;
}

std::string Like_KingPSFWithBkg::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: { std::ostringstream s; s << "r" << m_P*100; return s.str(); }
    case 1: return "1/gamma";
    case 2: return "P_s";
    case 3: return "N";
    default: throw std::string("Like_KingPSFWithBkg: unknown parameter");
    }
}

double Like_KingPSFWithBkg::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.25;
    case 1: return 0.5;
    case 2: return 0.5;
    case 3: return 1.0;
    default: throw std::string("Like_KingPSFWithBkg: unknown parameter");
    }  
}

double Like_KingPSFWithBkg::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.0;
    case 1: return 0.0;
    case 2: return 0.0;
    case 3: return 0.5;
    default: throw std::string("Like_KingPSFWithBkg: unknown parameter");
    }  
}

double Like_KingPSFWithBkg::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 10.0;
    case 1: return 1.0;
    case 2: return 2.0;
    case 3: return 2.0;
    default: throw std::string("Like_KingPSFWithBkg: unknown parameter");
    }  
}

double Like_KingPSFWithBkg::
modelRP(const std::vector<double>& p, const double P)
{
  const double rPm   = p[0];
  if(P==m_P)return rPm;
  const double gamma = 1.0/p[1];
  return 
    rPm*(pow(1.0-P,-1.0/(gamma-1.0))-1.0)/(pow(1.0-m_P,-1.0/(gamma-1.0))-1.0);
}

void Like_KingPSFWithBkg::
modelRP(double& rp, double& rp_err,
       const std::vector<double>& p,
       const std::vector<double>& p_err,
       const std::vector<std::vector<double> > p_cov,
       const double P)
{
  const double rPm   = p[0];
  if(P==m_P)
    {
      rp = rPm;
      rp_err = p_err[0];
      return;
    }

  rp = modelRP(p,P);

  const double drp_dp0 = rp/rPm;

  std::vector<double> pt(p);
  double dx = 0;
  dx = p[1]*1e-8;
  pt[1] += dx;
  const double drp_dp1 = (modelRP(pt,P)-rp)/dx;

  rp_err = std::sqrt(p_cov[0][0]*drp_dp0*drp_dp0
		     + p_cov[1][1]*drp_dp1*drp_dp1
		     + 2.0*p_cov[0][1]*drp_dp0*drp_dp1);
}

// ----------------------------------------------------------------------------
// 
// DOUBLE KING PROFILE
// 
// ----------------------------------------------------------------------------

Like_DblKingPSFWithBkg::
Like_DblKingPSFWithBkg(const std::vector<double>& d, double dmax,
			double P)
  : PSFCalc(), m_P(P), m_d2(), m_dmax2(dmax*dmax)
{
  const unsigned nd = d.size();
  m_d2.resize(nd);
  for(unsigned id=0;id<nd;id++)m_d2[id] = d[id]*d[id];
}

Like_DblKingPSFWithBkg::~Like_DblKingPSFWithBkg()
{
  // nothing to see here
}

unsigned Like_DblKingPSFWithBkg::numDim()
{
  return 1;
}

double Like_DblKingPSFWithBkg::Pdiff(const std::vector<double>& x,
					const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double rPm   = p[0];
  const double gamma1 = 1.0/p[1];
  const double gamma2 = 1.0/p[2];
  const double s2s1   = p[3];
  const double Ps2Ps  = p[5];
  const double sigma1 = _sigma(rPm,gamma1,gamma2,s2s1,Ps2Ps);
  const double sigma2 = s2s1*sigma1;

  const double Ps     = p[4];
  const double Ps1    = Ps*(1-Ps2Ps);
  const double Ps2    = Ps*Ps2Ps;
  const double Pb     = 1.0-Ps;
  const double N      = double(nd)*p[6];

  const double norm_s1 = N*Ps1*(gamma1-1)/(2*M_PI*sqr(sigma1))
    /(1-std::pow(1+m_dmax2/(2*sqr(sigma1)),1-gamma1));
  const double norm_s2 = N*Ps2*(gamma2-1)/(2*M_PI*sqr(sigma2))
    /(1-std::pow(1+m_dmax2/(2*sqr(sigma2)),1-gamma2));
  const double norm_b = N*Pb/(M_PI*m_dmax2);
  
  return norm_s1*std::pow(1+sqr(x[0])/(2*sqr(sigma1)),-gamma1) 
    + norm_s2*std::pow(1+sqr(x[0])/(2*sqr(sigma2)),-gamma2) 
    + norm_b;
}

double Like_DblKingPSFWithBkg::Pint(const std::vector<double>& x,
				       const std::vector<double>& x0,
				       const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double rPm   = p[0];
  const double gamma1 = 1.0/p[1];
  const double gamma2 = 1.0/p[2];
  const double s2s1   = p[3];
  const double Ps2Ps  = p[5];
  const double sigma1 = _sigma(rPm,gamma1,gamma2,s2s1,Ps2Ps);
  const double sigma2 = s2s1*sigma1;

  const double Ps     = p[4];
  const double Ps1    = Ps*(1-Ps2Ps);
  const double Ps2    = Ps*Ps2Ps;
  const double Pb     = 1.0-Ps;
  const double N      = double(nd)*p[6];

  const double norm_s1 = N*Ps1/(1-std::pow(1+m_dmax2/(2*sqr(sigma1)),1-gamma1));
  const double norm_s2 = N*Ps2/(1-std::pow(1+m_dmax2/(2*sqr(sigma2)),1-gamma2));
  const double norm_b = N*Pb/m_dmax2;

  return norm_s1*(std::pow(1+sqr(x0[0])/(2*sqr(sigma1)),1-gamma1) 
		  - std::pow(1+sqr(x[0])/(2*sqr(sigma1)),1-gamma1))
    + norm_s2*(std::pow(1+sqr(x0[0])/(2*sqr(sigma2)),1-gamma2) 
	       - std::pow(1+sqr(x[0])/(2*sqr(sigma2)),1-gamma2))
    + norm_b*(x[0]*x[0] - x0[0]*x0[0]);
}

std::string Like_DblKingPSFWithBkg::dimName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return "d";
    default: throw std::string("Like_DblKingPSFWithBkg: unknown dimension");
    }
}

bool Like_DblKingPSFWithBkg::
canCalcDFDP()
{
  return false;
}

double Like_DblKingPSFWithBkg::
f(const std::vector<double>& p)
{
  const unsigned nd = m_d2.size();

  const double rPm   = p[0];
  const double gamma1 = 1.0/p[1];
  const double gamma2 = 1.0/p[2];
  const double s2s1   = p[3];
  const double Ps2Ps  = p[5];
  const double sigma1 = _sigma(rPm,gamma1,gamma2,s2s1,Ps2Ps);
  const double sigma2 = s2s1*sigma1;

  const double Ps     = p[4];
  const double Ps1    = Ps*(1-Ps2Ps);
  const double Ps2    = Ps*Ps2Ps;
  const double Pb     = 1.0-Ps;
  const double N      = double(nd)*p[6];

  const double norm_s1 = N*Ps1*(gamma1-1)/(2*M_PI*sqr(sigma1))
    /(1-std::pow(1+m_dmax2/(2*sqr(sigma1)),1-gamma1));
  const double norm_s2 = N*Ps2*(gamma2-1)/(2*M_PI*sqr(sigma2))
    /(1-std::pow(1+m_dmax2/(2*sqr(sigma2)),1-gamma2));
  const double norm_b = N*Pb/(M_PI*m_dmax2);

  Accumulator loglike;
  for(unsigned id=0;id<nd;id++)
    loglike.add(std::log(norm_s1*std::pow(1+m_d2[id]/(2*sqr(sigma1)),-gamma1) 
			 + norm_s2*std::pow(1+m_d2[id]/(2*sqr(sigma2)),-gamma2) 
			 + norm_b));
  loglike.add(-N);
  return -loglike.sum();
}

double Like_DblKingPSFWithBkg::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("Like_DblKingPSFWithBkg: dfdp called");
}

unsigned Like_DblKingPSFWithBkg::numParam()
{
  return 7;
}

std::string Like_DblKingPSFWithBkg::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0: { std::ostringstream s; s << "r" << m_P*100; return s.str(); }
    case 1: return "1/gamma1";
    case 2: return "1/gamma2";
    case 3: return "sig2/sig1";
    case 4: return "P_s";
    case 5: return "Ps2/Ps";
    case 6: return "N";
    default: throw std::string("Like_DblKingPSFWithBkg: unknown parameter");
    }
}

double Like_DblKingPSFWithBkg::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.25;
    case 1: return 0.5;
    case 2: return 0.75;
    case 3: return 2.0;
    case 4: return 0.5;
    case 5: return 0.0;
    case 6: return 1.0;
    default: throw std::string("Like_DblKingPSFWithBkg: unknown parameter");
    }  
}

double Like_DblKingPSFWithBkg::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 0.0;
    case 1: return 0.0;
    case 2: return 0.0;
    case 3: return 0.0;
    case 4: return 0.0;
    case 5: return 0.0;
    case 6: return 0.5;
    default: throw std::string("Like_DblKingPSFWithBkg: unknown parameter");
    }  
}

double Like_DblKingPSFWithBkg::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0: return 10.0;
    case 1: return 1.0;
    case 2: return 1.0;
    case 3: return 1000.0;
    case 4: return 1.0;
    case 5: return 1.0;
    case 6: return 2.0;
    default: throw std::string("Like_DblKingPSFWithBkg: unknown parameter");
    }  
}

double Like_DblKingPSFWithBkg::
modelRP(const std::vector<double>& p, const double P)
{
  assert(0);
}

void Like_DblKingPSFWithBkg::
modelRP(double& rp, double& rp_err,
       const std::vector<double>& p,
       const std::vector<double>& p_err,
       const std::vector<std::vector<double> > p_cov,
       const double P)
{
  const double rPm   = p[0];
  if(P==m_P)
    {
      rp = rPm;
      rp_err = p_err[0];
      return;
    }
  assert(0);
}

double  Like_DblKingPSFWithBkg::
_sigma(double rP,double gamma1,double gamma2,double s2s1,double Ps2Ps)
{
  const double Ps1 = 1-Ps2Ps;
  const double Ps2 = Ps2Ps;
  const double rP2 = sqr(rP);

  //  std::cout << rP/sqrt(2.0*(pow(1.0-m_P,-1.0/(gamma1-1.0))-1.0)) << ' ' << rP/sqrt(2.0*(pow(1.0-m_P,-1.0/(gamma2-1.0))-1.0)) << '\n';

  double P = 1;
  double sr = 10;
  double sl = 0;
  while(P>m_P)
    {
      const double sigma1 = sr;
      const double sigma2 = sigma1*s2s1;
      P = Ps1*(1-std::pow(1+rP2/(2*sqr(sigma1)),-gamma1))
	+ Ps2*(1-std::pow(1+rP2/(2*sqr(sigma2)),-gamma2));
      if(P>m_P)sl = sr, sr *= 10;
    }
    
  unsigned nitteration = 0;
  P = 0;
  while(std::abs(P-m_P)/(1-m_P)>1e-6)
    {
      const double sm = 0.5*(sl+sr);
      const double sigma1 = sm;
      const double sigma2 = sigma1*s2s1;
      P = Ps1*(1-std::pow(1+rP2/(2*sqr(sigma1)),-gamma1))
	+ Ps2*(1-std::pow(1+rP2/(2*sqr(sigma2)),-gamma2));
      //std::cout << sl << ' ' << sr << ' ' << sm << ' ' << P << '\n';
      if(P>m_P)sl = sm;
      else sr = sm;
      nitteration++;
      assert(nitteration<1000);
    }

  return 0.5*(sl+sr);
}
