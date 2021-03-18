/** 
 * @file SplineFunction.h
 * @brief Declaration for the SplineFunction Function class
 * @author J. Chiang
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SplineFunction.h,v 1.3 2009/10/19 19:15:04 jchiang Exp $
 */

#ifndef Likelihood_SplineFunction_h
#define Likelihood_SplineFunction_h

#include "optimizers/Arg.h"
#include "optimizers/Function.h"

namespace Likelihood {

class RegularSpline;

/** 
 * @class SplineFunction
 *
 * @brief A function that uses a spline of N>=2 points between upper and lower
 * energy bounds. The spline is in log space [log10(E), log10(F)].
 *
 * @author S. Fegan
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/Likelihood/Likelihood/SplineFunction.h,v 1.3 2009/10/19 19:15:04 jchiang Exp $
 */

class SplineFunction : public optimizers::Function {

public:

   SplineFunction(double LowerLimit=100, double UpperLimit=100000);

   virtual double value(optimizers::Arg & x) const;

   virtual double derivByParam(optimizers::Arg & x, 
			       const std::string & paramName) const;

   // double integral(optimizers::Arg & xmin, optimizers::Arg & xmax) const;

   virtual void setAttributes(const Attributes& attributes);
   virtual void getAttributes(Attributes& attributes) const;
   virtual void initializeParameters();

   virtual Function * clone() const;

   unsigned npoints() const { return m_npoints; }

private:
   enum ParamTypes { LowerLimit, UpperLimit, Normalization, 
		     Index1, IndexN, Log10Flux0 };
   enum ControlPointPosition { CPP_EDGE, CPP_CENTER };
  
   SplineFunction(const SplineFunction& o);
   SplineFunction& operator=(const SplineFunction& o);

   void updateCache() const;

   bool                        m_initialized;
   unsigned                    m_npoints;
   ControlPointPosition        m_cpp;

   // cached variables
   mutable std::vector<double> m_cValues;
   mutable RegularSpline*      m_spline;
};

} // namespace Likelihood

#endif // Likelihood_SplineFunction_h
