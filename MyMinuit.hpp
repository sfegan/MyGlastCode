//-*-mode:c++; mode:font-lock;-*-

#ifndef MYMINUIT_HPP
#define MYMINUIT_HPP

#include <vector>
#include <string>
#include <map>

#include "Minuit/minuit_routines.h"

class Parameter
{
public:
  virtual ~Parameter();
  virtual std::string name() = 0;
};

class OptimizableParameter: public Parameter
{
  virtual ~OptimizableParameter();
  virtual double value() = 0;
  virtual void setValue(double x) = 0;
  virtual bool hasBounds() = 0;
  virtual double getLowBound() = 0;
  virtual double getHighBound() = 0;
  virtual double getDefaultValue() = 0;
};

class ParameterizableFunction
{
public:
  virtual ~ParameterizableFunction();
  virtual unsigned numParam() = 0;
  virtual std::string pName(unsigned iparam);
};

class ProbDensity: virtual public ParameterizableFunction
{
public:
  ProbDensity(): ParameterizableFunction() { }
  virtual ~ProbDensity();
  virtual unsigned numDim() = 0;
  virtual double Pdiff(const std::vector<double>& x,
		       const std::vector<double>& p) = 0;
};

class IntProbDensity: virtual ProbDensity
{
public:
  IntProbDensity(): ProbDensity() { }
  virtual ~IntProbDensity();
  virtual double Pint(const std::vector<double>& x,
		      const std::vector<double>& x0,
		      const std::vector<double>& p) = 0;
};

class Optimizable: virtual public ParameterizableFunction
{
public:
  Optimizable(): ParameterizableFunction() { }
  virtual ~Optimizable();
  virtual bool canCalcDFDP() = 0;
  virtual double f(const std::vector<double>& p) = 0;
  virtual double dfdp(const std::vector<double>& p, unsigned iparam) = 0;
  virtual double p0(unsigned iparam) = 0;
  virtual double plo(unsigned iparam) = 0;
  virtual double phi(unsigned iparam) = 0;
};

class ParameterModifier: public Optimizable
{
public:
  ParameterModifier(Optimizable* fn);
  virtual ~ParameterModifier();
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
  void fixPVal(unsigned iparam, double pval);
  void setP0(unsigned iparam, double pval);
  void setPLo(unsigned iparam, double pval);
  void setPHi(unsigned iparam, double pval);
private:
  Optimizable* m_fn;
  std::map<std::string, unsigned> m_pmap;
  std::map<unsigned, double> m_pfixed;
  std::map<unsigned, double> m_p0;
  std::map<unsigned, double> m_phi;
  std::map<unsigned, double> m_plo;
  unsigned iFnParam(unsigned iparam) const;
};

class MyMinuit
{
public:
  MyMinuit(Optimizable* fn, int verbose = 1); 
  void setEpsilon(double epsilon);
  void minimize(bool interactive = false, double tol=1e-3, double err=0.5,
		const std::string& command_str = "", bool nohesse=false,
		unsigned strategy = 1);

  // Parameter names and numbers
  unsigned numParam() const { return m_fn->numParam(); }
  std::vector<std::string> pNames() const;
  std::string pName(unsigned iparam) const;

  /// Function and parameter values
  double fVal() const;
  std::vector<double> pVals() const;
  double pVal(unsigned iparam) const;
  std::vector<double> pErrs() const;
  double pErr(unsigned iparam) const;
  std::vector<std::vector<double> > pCovMtx() const;

  struct MinStat
  {
    double fmin;    //! Best function value found so far
    double fedm;    //! Estimated vertical distance remaining to minimum
    double errdef;  //! Value of UP, defining parameter uncertainties
    int npari;      //! Number of free parameters
    int nparx;      //! Number of total parameters
    int istat;      //! Status indicating quality of covariance matrix
  };
    
  /// Get the current status of the minimization
  MinStat getMinimizationStatus() const;

  struct ValAndErr
  {
    std::string pname;    //! Short name
    std::string mname;    //! Minuit name
    double val;           //! Value
    double err;           //! Error estimate 
    double bnd_lo;        //! Left boundary
    double bnd_hi;        //! Right boundary
    double err_plus;      //! Error from MINOS (plus)
    double err_minus;     //! Error from MINOS (minus)
    double err_parabolic; //! Parabolic error estimate
    double global_corr;   //! Some weird thing
  };

  ValAndErr getValueAndError(unsigned iparam) const;
  std::vector<ValAndErr> getValuesAndErrors() const;

  int computeMinosErrors(unsigned iparam, 
			 double& perr_minus, double& perr_plus);
  int computeMinosContour(unsigned xparam, unsigned yparam,
			  std::vector<double>& x, std::vector<double>& y,
			  unsigned n = 20);

private:
  void minimization_program_fixed(double tol, bool nohesse, unsigned strategy);
  void minimization_program_interactive(double tol, 
					const std::string& command_str);

  int doCmd(const std::string& command);

  void fcn(int* npar, double* grad, double* fcnval,
	   double* pval, int* iflag);

  static void s_fcn(int* npar, double* grad, double* fcnval,
		    double* pval, int* iflag, void* futil);

  //! Split command text into words and perform variable subsitution
  std::vector<std::string> 
  splitCommand(const std::string& cmdtxt, 
	       const std::map<std::string,unsigned>& parameters,
	       double tol) const;

  //! Print help for interactive commands
  void printHelp(bool long_help) const;

  //! Print parameters
  void printParameters() const;

  //! Print dump of current state
  void dumpState(std::ostream& str, bool dumpCovarianceMatrix=true) const;

  //! Functions relating to chi^2 calculation
  static double invchisquaredistribution(double v, double y);
  static double invincompletegammac(double a, double y0);
  static double incompletegammac(double a, double x);
  static double incompletegamma(double a, double x);
  static double lngamma(double x, double& sgngam);
  static double invnormaldistribution(double y0);
  
  Optimizable* m_fn;
  int m_verbose;
};

#endif
