//-*-mode:c++; mode:font-lock;-*-

/*! \file VSASVD.hpp

  Singular Value Decomposition

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/02/2007
*/

#ifndef VSASVD_HPP
#define VSASVD_HPP

#include <limits>

#include <VSAAlgebra.hpp>

#define VSA_ASSERT

namespace VERITAS
{
  namespace VSAAlgebra
  {

    class SVD
    {
    public:
      SVD(const VSAAlgebra::MatrixND& a);
    
      const VSAAlgebra::MatrixND& u() const { return m_u; }
      const VSAAlgebra::MatrixND& v() const { return m_v; }
      const VSAAlgebra::VecND& w() const { return m_w; }
      double thresh() const { return m_thresh; }

      void solve(const VSAAlgebra::VecND& b, VSAAlgebra::VecND& x, 
		 double thresh = -1.0);
      void solve(const VSAAlgebra::MatrixND& b, VSAAlgebra::MatrixND& x, 
		 double thresh = -1.0);

      unsigned rank(double thresh = -1);
      unsigned nullity(double thresh = -1);

      void range(VSAAlgebra::MatrixND& m, double thresh = -1);
      void nullSpace(VSAAlgebra::MatrixND& m, double thresh = -1);

      double invCondition() 
      { return (m_w[0]<=0 || m_w[m_ncol-1]<=0)?0.0:m_w[m_ncol-1]/m_w[0]; }

    private:
      void decompose();
      void reorder();
      double pythag(const double a, const double b);
      void setThresh(double thresh);

      int                    m_nrow;
      int                    m_ncol;
      VSAAlgebra::MatrixND   m_u;
      VSAAlgebra::MatrixND   m_v;
      VSAAlgebra::VecND      m_w;

      double                 m_eps;
      double                 m_thresh;
    };

  }
}

#endif // ifndef VSASVD_HPP
