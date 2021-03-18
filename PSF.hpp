//-*-mode:c++; mode:font-lock;-*-

#ifndef PSF_HPP
#define PSF_HPP

#include <vector>

#include "MyMinuit.hpp"

class PSFCalc
  : virtual public Optimizable, virtual public IntProbDensity
{
public:
  PSFCalc(): Optimizable(), IntProbDensity() { }
  virtual ~PSFCalc();
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68) = 0;
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68) = 0;
};

class Like_IsotropicBkg: public PSFCalc
{
public:
  Like_IsotropicBkg(const std::vector<double>& d, double dmax);
  virtual ~Like_IsotropicBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  double m_nd;
  double m_dmax2;
};

class Like_GaussianPSFWithBkg: public PSFCalc
{
public:
  Like_GaussianPSFWithBkg(const std::vector<double>& d, double dmax);
  virtual ~Like_GaussianPSFWithBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  std::vector<double> m_d2;
  double m_dmax2;
};

class Like_DblGaussianPSFWithBkg: public PSFCalc
{
public:
  Like_DblGaussianPSFWithBkg(const std::vector<double>& d, double dmax);
  virtual ~Like_DblGaussianPSFWithBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  std::vector<double> m_d2;
  double m_dmax2;
};

class Like_MultiGaussianPSFWithBkg: public PSFCalc
{
public:
  Like_MultiGaussianPSFWithBkg(unsigned n,
			       const std::vector<double>& d, double dmax);
  virtual ~Like_MultiGaussianPSFWithBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  double I(double x, const std::vector<double>& p) const;
  unsigned m_n;
  std::vector<double> m_d2;
  double m_dmax2;
};

class Like_BurdettPSFWithBkg: public PSFCalc
{
public:
  Like_BurdettPSFWithBkg(const std::vector<double>& d, double dmax);
  virtual ~Like_BurdettPSFWithBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  std::vector<double> m_d2;
  double m_dmax2;
};

class Like_KingPSFWithBkg: public PSFCalc
{
public:
  Like_KingPSFWithBkg(const std::vector<double>& d, double dmax, 
		       double P=0.68);
  virtual ~Like_KingPSFWithBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  double _sigma(double rP, double gamma) 
  {
    return rP/sqrt(2.0*(pow(1.0-m_P,-1.0/(gamma-1.0))-1.0));
  }

  double m_P;
  std::vector<double> m_d2;
  double m_dmax2;
};

class Like_DblKingPSFWithBkg: public PSFCalc
{
public:
  Like_DblKingPSFWithBkg(const std::vector<double>& d, double dmax, 
			  double P=0.68);
  virtual ~Like_DblKingPSFWithBkg();
  virtual unsigned numDim();
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p);
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p);
  virtual std::string dimName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  virtual double modelRP(const std::vector<double>& p,
			const double P = 0.68);
  virtual void modelRP(double& rp, double& rp_err,
		       const std::vector<double>& p,
		       const std::vector<double>& p_err,
		       const std::vector<std::vector<double> > p_cov,
		       const double P = 0.68);
private:
  double _sigma(double rP,double gamma1,double gamma2,double s2s1,double Ps2Ps);

  double m_P;
  std::vector<double> m_d2;
  double m_dmax2;
};

#endif
