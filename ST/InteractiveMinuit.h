/** 
 * @file InteractiveMinuit.h
 * @brief Declaration for the Interactive Minuit Optimizer subclass.
 * @author P. Nolan
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/optimizers/Minuit.h,v 1.14 2008/09/29 16:08:17 jchiang Exp $
 */

#ifndef optimizers_INTERACTIVE_MINUIT_H
#define optimizers_INTERACTIVE_MINUIT_H

#include <map>
#include "optimizers/Optimizer.h"
#include "optimizers/Statistic.h"
#include "optimizers/f2c_types.h"

namespace optimizers {

  /** 
   * @class InteractiveMinuit
   *
   * @brief Wrapper class for the Minuit optimizer from CERN
   *
   * @author P. Nolan
   *    
   This class implements an Optimizer by using Minuit, a well-known
   package from CERN.  All variables are treated as bounded.  User
   interaction is required.
   */
  
  // Doxygen the C file here so it can be left as nearly as
  // possible in its pristine, machine-produced state.
  /**
   * @file minuit_routines.c
   *
   * @brief The Minuit package translated from Fortran by f2c
   *
   Minuit is a well-known optimizing/fitting package in the HEP
   community.  It has been developed for over 30 years at CERN.
   
   This file was produced from the CERN Fortran source code.
   First, g77 -E was used to insert all the \#include files and
   make a single, large Fortran file free of preprocessor commands.
   Then f2c -C++ produced this file.  The only hand modification
   required was to change \#include "f2c.h" to \#include "f2c/f2c.h"
   to conform to the way CMT wants files to be laid out.
   
   In non-interactive mode, the API for using Minuit involves the
   functions mninit_, mnseti_, mnparm_, mnpars_, mnexcm_, mncomd_,
   mnpout_, mnstat_, mnemat_, mnerrs_, mncont_, mnintr_, and mninpu_.
  */
  
  class InteractiveMinuit : public Optimizer {
    
  public:
    
    InteractiveMinuit(Statistic &stat, const std::string& cmdtxt);
    
    virtual ~InteractiveMinuit() {}
    
    virtual int find_min(int verbose, double tol, int tolType = ABSOLUTE);
    virtual int find_min_only(int verbose, double tol, int tolType = ABSOLUTE);

    //! Run a MINOS error analysis
    std::pair<double,double> Minos(unsigned int n);

    //! One-sigma confidence regions based on Hessian, assuming 
    // that this function is a likelihood
    virtual const std::vector<double> & getUncertainty(bool useBase=false);

    /// Access to the covariance matrix
    virtual std::vector< std::vector<double> > covarianceMatrix() const;

    /// Doesn't do very much here
    virtual std::ostream& put(std::ostream& s) const;

    /// Get list of parameter names
    std::vector<std::string> getParameterList() const;
    
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
      std::string uname;    //! Unique name
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

    std::vector<ValAndErr> getValuesAndErrors() const;

    std::vector<std::vector<double> > getCovarianceMatrix() const;

    int computeMinosContour(unsigned xparam, unsigned yparam,
                            std::vector<double>& x, std::vector<double>& y,
                            unsigned n = 20);

  private:
    
    //! Pass a command string to Minuit
    int doCmd(std::string command);
    
    //! The function which InteractiveMinuit will minimize
    static void fcn(int* npar, double* grad, double* fcnval,
		    double* xval, int* iflag, void* futil);

    //! Split command text into words and perform variable subsitution
    std::vector<std::string> 
      splitCommand(const std::string& cmdtxt, 
		   const std::map<std::string,unsigned>& parameters,
		   double tolerance) const;

    //! Print help for interactive commands
    void printHelp(bool long_help) const;

    //! Print parameters
    void printParameters() const;

    //! Print dump of current state
    void dumpState(std::ostream& str, bool dumpCovarianceMatrix=true) const;
    void unDumpState(std::istream& str);

    //! Functions relating to chi^2 calculation
    static double invchisquaredistribution(double v, double y);
    static double invincompletegammac(double a, double y0);
    static double incompletegammac(double a, double x);
    static double incompletegamma(double a, double x);
    static double lngamma(double x, double& sgngam);
    static double invnormaldistribution(double y0);

    std::vector<std::string> m_commands;
    int                      m_nparam;
    int                      m_quality;
    double                   m_val;
    double                   m_distance;
    std::vector<double>      m_stored_param_values;
  };
  
} // namespace optimizers


#endif // optimizers_INTERACTIVE_MINUIT_H
