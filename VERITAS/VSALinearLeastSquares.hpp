//-*-mode:c++; mode:font-lock;-*-

/*! \file VSALinearLeastSquares.hpp

  Generalized Linear Least Squares Fitting

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2007
*/

#ifndef VSALINEARLEASTSQUARES_HPP
#define VSALINEARLEASTSQUARES_HPP

#include <vector>
#include <stdexcept>

#include <VSAAlgebra.hpp>
#include <VSASVD.hpp>

#define VSA_ASSERT

namespace VERITAS
{

  namespace VSAMath
  {

    struct DataPoint
    {
      DataPoint(): x(), y(), sigma() { /* nothing to see here */ }
      DataPoint(double _x, double _y, double _s=1.0): 
	x(_x), y(_y), sigma(_s) { /* nothing to see here */ }

      double x;
      double y;
      double sigma;
    };

    typedef std::vector<DataPoint> Data;

    template<typename Fn> class Fitlin
    {
    public:
      struct Options { /* nothing to see here */ };
      static Options defaultOptions() { return Options(); }

      Fitlin(const Data& data, Fn& fn, const Options& opt = defaultOptions()):
	m_ndata(data.size()), m_data(data), m_fn(fn), 
	m_nparm(), m_fit(), m_a(), m_cov(), m_chi2()
      {
	VSAAlgebra::VecND atemp;
	m_fn(m_data[0].x, atemp);

	m_nparm = atemp.ndim();
	m_a     .resize(m_nparm); 
	m_fit   .resize(m_nparm,true); 
	m_cov   .resize(m_nparm,m_nparm); 
      }


      void hold(unsigned i, double val) { m_fit[i] = false; m_a[i]=val; }
      void free(unsigned i) { m_fit[i] = true; }
      
      void fit();

      const VSAAlgebra::VecND& param() const { return m_a; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }

    private:
      unsigned               m_ndata;
      Data                   m_data;
      Fn                     m_fn;

      unsigned               m_nparm;
      std::vector<bool>      m_fit;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::MatrixND   m_cov;
      double                 m_chi2;
    };

    template<typename Fn> class Fitsvd
    {
    public:
      typedef double Options;
      static Options defaultOptions() { return 1e-12; }

      Fitsvd(const Data& data, Fn& fn, const double& tol = defaultOptions()):
	m_ndata(data.size()), m_data(data), m_fn(fn), m_tol(tol),
	m_nparm(), m_a(), m_cov(), m_chi2()
      {
	VSAAlgebra::VecND atemp;
	m_fn(m_data[0].x, atemp);

	m_nparm = atemp.ndim();
	m_a     .resize(m_nparm); 
	m_cov   .resize(m_nparm,m_nparm); 
      }

      double tol() const { return m_tol; }
      void setTol(const double& tol) { m_tol = tol; }

      void fit();

      const VSAAlgebra::VecND& param() const { return m_a; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }

    private:
      unsigned               m_ndata;
      Data                   m_data;
      Fn                     m_fn;
      double                 m_tol;

      unsigned               m_nparm;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::MatrixND   m_cov;
      double                 m_chi2;
    };

    template<typename Fn> void Fitlin<Fn>::fit()
    {
      // Based on NR3 algorithm - Section 15.4
      
      unsigned nfit = 0;
      for(unsigned iparm=0;iparm<m_nparm;iparm++)if(m_fit[iparm])nfit++;
      if(nfit==0)throw
	std::domain_error(std::string(__PRETTY_FUNCTION__)+
			  ": no parameters to be fitted");
      
      VSAAlgebra::MatrixND temp(nfit,nfit,0.0);
      VSAAlgebra::MatrixND beta(nfit,1,0.0);

      VSAAlgebra::VecND yfunc(m_nparm);      
      for(unsigned idata=0;idata<m_ndata;idata++)
	{
	  m_fn(m_data[idata].x, yfunc);
	  double ym = m_data[idata].y;
	  double sm = m_data[idata].sigma;
	  if(sm<=0)throw 
	    std::domain_error(std::string(__PRETTY_FUNCTION__) +
			      ": RMS of data point is zero or negative");

	  if(nfit < m_nparm)
	    for(unsigned iparm=0;iparm<m_nparm;iparm++)
	      if(!m_fit[iparm])ym -= m_a[iparm]*yfunc[iparm];

	  double sig2inv = 1.0/(sm*sm);
	  for(unsigned iparm=0, ifit=0; iparm<m_nparm; iparm++)
	    if(m_fit[iparm])
	      {
		const double wt = yfunc[iparm]*sig2inv;
		for(unsigned jparm=0, jfit=0; jparm<m_nparm; jparm++)
		  if(m_fit[jparm])temp(ifit,jfit++) += wt*yfunc[jparm];
		beta(ifit++,0) += wt*ym;
	      }
	}

      for(unsigned ifit=1; ifit<nfit; ifit++)
	for(unsigned jfit=0; jfit<ifit; jfit++)
	  temp(jfit,ifit) = temp(ifit,jfit);
      
      VSAAlgebra::MatrixND::gaussJordan(temp, beta);

      for(unsigned iparm=0, ifit=0; iparm<m_nparm; iparm++)
	if(m_fit[iparm])m_a[iparm] = beta(ifit++,0);

      m_chi2 = 0.0;
      for(unsigned idata=0; idata<m_ndata; idata++)
	{
	  m_fn(m_data[idata].x, yfunc);
	  double sum = 0.0;
	  for(unsigned iparm=0; iparm<m_nparm; iparm++)
	    sum += m_a[iparm]*yfunc[iparm];
	  sum = (m_data[idata].y - sum)/m_data[idata].sigma;
	  m_chi2 += sum*sum;
	}

      for(unsigned ifit=0; ifit<nfit; ifit++)
	for(unsigned jfit=0; jfit<nfit; jfit++)
	  m_cov(ifit, jfit) = temp(ifit, jfit);

      for(unsigned iparm=nfit; iparm<m_nparm; iparm++)
	for(unsigned jparm=0; jparm<iparm+1; jparm++)
	  m_cov(iparm,jparm) = m_cov(jparm,iparm) = 0.0;
      
      unsigned kparm = nfit-1;
      for(int jparm = m_nparm-1; jparm>=0; jparm--)
	if(m_fit[jparm])
	  {
	    for(unsigned iparm=0; iparm<m_nparm; iparm++)
	      std::swap(m_cov(iparm,kparm), m_cov(iparm, jparm));
	    for(unsigned iparm=0; iparm<m_nparm; iparm++)
	      std::swap(m_cov(kparm,iparm), m_cov(jparm, iparm));
	    kparm--;
	  }
    }

    template<typename Fn> void Fitsvd<Fn>::fit()
    {
      // Based on NR3 algorithm - Section 15.4
      
      VSAAlgebra::MatrixND aa(m_ndata,m_nparm);
      VSAAlgebra::VecND b(m_ndata);
      VSAAlgebra::VecND yfunc(m_nparm);      
      for(unsigned idata=0;idata<m_ndata;idata++)
	{
	  m_fn(m_data[idata].x, yfunc);
	  const double ym = m_data[idata].y;
	  const double sm = m_data[idata].sigma;
	  if(sm<=0)throw 
	    std::domain_error(std::string(__PRETTY_FUNCTION__) +
			      ": RMS of data point is zero or negative");

	  const double siginv = 1.0/sm;
	  for(unsigned iparm=0; iparm<m_nparm; iparm++)
	    aa(idata,iparm) = yfunc(iparm)*siginv;
	  b[idata] = ym*siginv;
	}

      VSAAlgebra::SVD svd(aa);
      double thresh = (m_tol > 0 ? m_tol*svd.w()[0] : -1);
      svd.solve(b,m_a,thresh);

      m_chi2 = 0.0;
      for(unsigned idata=0; idata<m_ndata; idata++)
	{
	  m_fn(m_data[idata].x, yfunc);
	  double sum = 0.0;
	  for(unsigned iparm=0; iparm<m_nparm; iparm++)
	    sum += aa(idata,iparm)*m_a[iparm];
	  sum -= b[idata];
	  m_chi2 += sum*sum;
	}

      for(unsigned iparm=0; iparm<m_nparm; iparm++)
	for(unsigned jparm=0; jparm<iparm+1; jparm++)
	  {
	    double sum = 0;
	    for(unsigned kparm=0; kparm<m_nparm; kparm++)
	      {
		const double wk = svd.w()[kparm];
		if(wk > svd.thresh())
		  sum += svd.v()(iparm,kparm)*svd.v()(jparm,kparm)/(wk*wk);
	      }
	    m_cov(jparm,iparm) = m_cov(iparm,jparm) = sum;
	  }
    }

    class PolyFit
    {
    public:

      // ----------------------------------------------------------------------
      // Functor used by fitter to calculate value of basis functions at data
      // x-values. Returns value of th n+1 polynomial functions (x^0 .. x^n)
      // in a VSAAlgebra::VecND
      // ----------------------------------------------------------------------

      class Fn
      {
      public:
	Fn(unsigned n): m_n(n+1) { /* nothing to see here */ }
	void operator() (double x, VSAAlgebra::VecND& v)
	{
	  v.resize(m_n);
	  double s = v[0] = 1.0;
	  for(unsigned i=1;i<m_n;i++)v[i] = s *= x;
	}
      private:
	unsigned m_n;
      };

      // ----------------------------------------------------------------------
      // Static interface class to free user from need to manually construct 
      // the fitter, call it and extract the results. Simply call
      // "PolyFit::fit" and let it do all the work
      // ----------------------------------------------------------------------

      typedef Fitsvd<Fn> PreferredFitter;

      static double fit(unsigned n, const Data& data, VSAAlgebra::VecND& param,
			VSAAlgebra::MatrixND* cov = 0, 
			const PreferredFitter::Options& opt =
			PreferredFitter::defaultOptions());

      static VSAAlgebra::VecND fit(unsigned n, const Data& data)
      {
	VSAAlgebra::VecND param;
	fit(n, data, param);
	return param;
      }

      static double val(const VSAAlgebra::VecND& a, double x);

      static void val(const VSAAlgebra::VecND& a, const std::vector<double> x,
		      std::vector<double>& y);

      static std::vector<double> 
      val(const VSAAlgebra::VecND& a, const std::vector<double> x)
      {
	std::vector<double> y;
	val(a,x,y);
	return y;
      }

      static void differentiate(const VSAAlgebra::VecND& a, 
				VSAAlgebra::VecND& dydx_a);

      static VSAAlgebra::VecND differentiate(const VSAAlgebra::VecND& a)
      {
	VSAAlgebra::VecND dydx_a;
	differentiate(a, dydx_a);
	return dydx_a;
      }

    };

  }
}

#endif // #ifndef VSALINEARLEASTSQUARES_HPP
