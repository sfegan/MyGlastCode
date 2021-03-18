//-*-mode:c++; mode:font-lock;-*-

/*! \file VSASVD.cpp

  Singular Value Decomposition

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    0.1
  \date       12/02/2007
*/

#include <stdexcept>
#include <VSASVD.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

SVD::SVD(const MatrixND& a):
  m_nrow(a.nrow()), m_ncol(a.ncol()), m_u(a), m_v(m_ncol,m_ncol), m_w(m_ncol),
  m_eps(std::numeric_limits<double>::epsilon()), m_thresh()
{
  decompose();
  reorder();
  setThresh(-1.0);
}

void SVD::solve(const VecND& b, VecND& x, double thresh)
{
  if(b.ndim() != unsigned(m_nrow) || x.ndim() != unsigned(m_ncol))throw
    std::invalid_argument(std::string(__PRETTY_FUNCTION__)
			  + ": bad input vector sizes");

  VecND temp(m_ncol); //m_nrow);
  setThresh(thresh);
  
  for(int jcol=0;jcol<m_ncol;jcol++)
    {
      double s = 0.0;
      if(m_w[jcol] > m_thresh)
	{
	  for(int irow=0;irow<m_nrow;irow++)s += m_u(irow,jcol)*b[irow];
	  s /= m_w[jcol];
	}
      temp[jcol] = s;
    }
  
  for(int jcol=0;jcol<m_ncol;jcol++)
    {
      double s = 0.0;
      for(int jjcol = 0;jjcol<m_ncol;jjcol++)
	s += m_v(jcol,jjcol)*temp[jjcol];
      x[jcol] = s;
    }
}

void SVD::solve(const MatrixND& b, MatrixND& x, double thresh)
{
  if(b.nrow() != unsigned(m_ncol) || x.nrow() != unsigned(m_ncol)
     || b.ncol() != x.ncol())
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__)
				+ ": bad input matrix sizes");
  
  VecND xx(m_ncol);
  for(int jrow=0; jrow<m_nrow; jrow++)
    {
      for(int icol=0; icol<m_ncol; icol++)xx[icol] = b(icol,jrow);
      solve(xx,xx,thresh);
      for(int icol=0; icol<m_ncol; icol++)x(icol,jrow) = xx[icol];
    }
}

unsigned SVD::rank(double thresh)
{
  setThresh(thresh);
  unsigned n = 0;
  for(int jcol=0;jcol<m_ncol;jcol++)if(m_w[jcol] > m_thresh)n++;
  return n;
}

unsigned SVD::nullity(double thresh)
{
  setThresh(thresh);
  unsigned n = 0;
  for(int jcol=0;jcol<m_ncol;jcol++)if(m_w[jcol] <= m_thresh)n++;
  return n;
}

void SVD::range(VSAAlgebra::MatrixND& m, double thresh)
{
  m.resize(m_nrow,rank(thresh));
  unsigned n=0;
  for(int jcol=0;jcol<m_ncol;jcol++)
    if(m_w[jcol] > m_thresh)
      {
	for(int irow=0;irow<m_nrow;irow++)
	  m(irow,n) = m_u(irow,jcol);
	n++;
      }
}

void SVD::nullSpace(VSAAlgebra::MatrixND& m, double thresh)
{
  m.resize(m_nrow,nullity(thresh));
  unsigned n=0;
  for(int jcol=0;jcol<m_ncol;jcol++)
    if(m_w[jcol] <= m_thresh)
      {
	for(int irow=0;irow<m_nrow;irow++)
	  m(irow,n) = m_u(irow,jcol);
	n++;
      }
}

void SVD::setThresh(double thresh)
{
  if(thresh>=0.0)m_thresh = thresh;
  else m_thresh = 0.5*sqrt(double(m_nrow+m_ncol)+1.0)*m_w[0]*m_eps;
}

// ============================================================================
//
// INTERNAL FUNCTIONS: Straight from NR3 Webnote 2 with some REGEXP changes
//
// ============================================================================

inline static double SIGN(const double& a, const double& b)
{
  return b >= 0.0 ? fabs(a) : -fabs(a);
}

void SVD::decompose() 
{
  bool flag;
  int i,its,j,jj,k,l=0,nm=0;
  double anorm,c,f,g,h,s,scale,x,y,z;
  VecND rv1(m_ncol);
  g = scale = anorm = 0.0;
  for (i=0;i<m_ncol;i++) 
    {
      l=i+2;
      rv1[i]=scale*g;
      g=s=scale=0.0;
      if (i < m_nrow) 
	{
	  for (k=i;k<m_nrow;k++) scale += fabs(m_u(k,i));
	  if (scale != 0.0) 
	    {
	      for (k=i;k<m_nrow;k++) 
		{
		  m_u(k,i) /= scale;
		  s += m_u(k,i)*m_u(k,i);
		}
	      f=m_u(i,i);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      m_u(i,i)=f-g;
	      for (j=l-1;j<m_ncol;j++) 
		{
		  for (s=0.0,k=i;k<m_nrow;k++) s += m_u(k,i)*m_u(k,j);
		  f=s/h;
		  for (k=i;k<m_nrow;k++) m_u(k,j) += f*m_u(k,i);
		}
	      for (k=i;k<m_nrow;k++) m_u(k,i) *= scale;
	    }
	}
      m_w[i]=scale *g;
      g=s=scale=0.0;
      if (i+1 <= m_nrow && i+1 != m_ncol) 
	{
	  for (k=l-1;k<m_ncol;k++) scale += fabs(m_u(i,k));
	  if (scale != 0.0) 
	    {
	      for (k=l-1;k<m_ncol;k++) 
		{
		  m_u(i,k) /= scale;
		  s += m_u(i,k)*m_u(i,k);
		}
	      f=m_u(i,l-1);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      m_u(i,l-1)=f-g;
	      for (k=l-1;k<m_ncol;k++) rv1[k]=m_u(i,k)/h;
	      for (j=l-1;j<m_nrow;j++) 
		{
		  for (s=0.0,k=l-1;k<m_ncol;k++) s += m_u(j,k)*m_u(i,k);
		  for (k=l-1;k<m_ncol;k++) m_u(j,k) += s*rv1[k];
		}
	      for (k=l-1;k<m_ncol;k++) m_u(i,k) *= scale;
	    }
	}
      anorm=std::max(anorm,(fabs(m_w[i])+fabs(rv1[i])));
    }

  for (i=m_ncol-1;i>=0;i--)
    { 
      if (i < m_ncol-1)
	{
	  if (g != 0.0) 
	    {
	      for (j=l;j<m_ncol;j++)
		m_v(j,i)=(m_u(i,j)/m_u(i,l))/g;
	      for (j=l;j<m_ncol;j++) 
		{
		  for (s=0.0,k=l;k<m_ncol;k++) s += m_u(i,k)*m_v(k,j);
		  for (k=l;k<m_ncol;k++) m_v(k,j) += s*m_v(k,i);
		}
	    }
	  for (j=l;j<m_ncol;j++) m_v(i,j)=m_v(j,i)=0.0;
	}
      m_v(i,i)=1.0;
      g=rv1[i];
      l=i;
    }

  for (i=std::min(m_nrow,m_ncol)-1;i>=0;i--) 
    {
      l=i+1;
      g=m_w[i];
      for (j=l;j<m_ncol;j++) m_u(i,j)=0.0;
      if (g != 0.0) 
	{
	  g=1.0/g;
	  for (j=l;j<m_ncol;j++) 
	    {
	      for (s=0.0,k=l;k<m_nrow;k++) s += m_u(k,i)*m_u(k,j);
	      f=(s/m_u(i,i))*g;
	      for (k=i;k<m_nrow;k++) m_u(k,j) += f*m_u(k,i);
	    }
	  for (j=i;j<m_nrow;j++) m_u(j,i) *= g;
	} 
      else for (j=i;j<m_nrow;j++) m_u(j,i)=0.0;
      ++m_u(i,i);
    }

  for (k=m_ncol-1;k>=0;k--) 
    {
      for (its=0;its<30;its++)
	{
	  flag=true;
	  for (l=k;l>=0;l--) 
	    {
	      nm=l-1;
	      if (l == 0 || fabs(rv1[l]) <= m_eps*anorm)
		{
		  flag=false;
		  break;
		}
	      if (fabs(m_w[nm]) <= m_eps*anorm) break;
	    }
	  if (flag) 
	    {
	      c=0.0;
	      s=1.0;
	      for (i=l;i<k+1;i++) 
		{
		  f=s*rv1[i];
		  rv1[i]=c*rv1[i];
		  if (fabs(f) <= m_eps*anorm) break;
		  g=m_w[i];
		  h=pythag(f,g);
		  m_w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s = -f*h;
		  for (j=0;j<m_nrow;j++) 
		    {
		      y=m_u(j,nm);
		      z=m_u(j,i);
		      m_u(j,nm)=y*c+z*s;
		      m_u(j,i)=z*c-y*s;
		    }
		}
	    }
	  z=m_w[k];
	  if (l == k) 
	    { 
	      if (z < 0.0) 
		{
		  m_w[k] = -z;
		  for (j=0;j<m_ncol;j++) m_v(j,k) = -m_v(j,k);
		}
	      break;
	    }
	  if (its == 29) 
	    throw std::overflow_error(std::string(__PRETTY_FUNCTION__)
				  + ": no convergence after 30 iterations");
	  x=m_w[l];
	  nm=k-1;
	  y=m_w[nm];
	  g=rv1[nm];
	  h=rv1[k];
	  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g=pythag(f,1.0);
	  f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c=s=1.0; 
	  for (j=l;j<=nm;j++) 
	    {
	      i=j+1;
	      g=rv1[i];
	      y=m_w[i];
	      h=s*g;
	      g=c*g;
	      z=pythag(f,h);
	      rv1[j]=z;
	      c=f/z;
	      s=h/z;
	      f=x*c+g*s;
	      g=g*c-x*s;
	      h=y*s;
	      y *= c;
	      for (jj=0;jj<m_ncol;jj++) 
		{
		  x=m_v(jj,j);
		  z=m_v(jj,i);
		  m_v(jj,j)=x*c+z*s;
		  m_v(jj,i)=z*c-x*s;
		}
	      z=pythag(f,h);
	      m_w[j]=z;
	      if (z) 
		{
		  z=1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f=c*g+s*y;
	      x=c*y-s*g;
	      for (jj=0;jj<m_nrow;jj++) 
		{
		  y=m_u(jj,j);
		  z=m_u(jj,i);
		  m_u(jj,j)=y*c+z*s;
		  m_u(jj,i)=z*c-y*s;
		}
	    }
	  rv1[l]=0.0;
	  rv1[k]=f;
	  m_w[k]=x;
	}
    }
}

void SVD::reorder() 
{
  int i,j,k,s,inc=1;
  double sw;
  VecND su(m_nrow);
  VecND sv(m_ncol);
  do { inc *= 3; inc++; } while (inc <= m_ncol);
  do {
    inc /= 3;
    for (i=inc;i<m_ncol;i++) {
      sw = m_w[i];
      for (k=0;k<m_nrow;k++) su[k] = m_u(k,i);
      for (k=0;k<m_ncol;k++) sv[k] = m_v(k,i);
      j = i;
      while (m_w[j-inc] < sw) {
	m_w[j] = m_w[j-inc];
	for (k=0;k<m_nrow;k++) m_u(k,j) = m_u(k,j-inc);
	for (k=0;k<m_ncol;k++) m_v(k,j) = m_v(k,j-inc);
	j -= inc;
	if (j < inc) break;
      }
      m_w[j] = sw;
      for (k=0;k<m_nrow;k++) m_u(k,j) = su[k];
      for (k=0;k<m_ncol;k++) m_v(k,j) = sv[k];
    }
  } while (inc > 1);
  for (k=0;k<m_ncol;k++) {
    s=0;
    for (i=0;i<m_nrow;i++) if (m_u(i,k) < 0.) s++;
    for (j=0;j<m_ncol;j++) if (m_v(j,k) < 0.) s++;
    if (s > (m_nrow+m_ncol)/2) {
      for (i=0;i<m_nrow;i++) m_u(i,k) = -m_u(i,k);
      for (j=0;j<m_ncol;j++) m_v(j,k) = -m_v(j,k);
    }
  }
}

double SVD::pythag(const double a, const double b) 
{
  double fabsa=fabs(a), fabsb=fabs(b);
  return (fabsa > fabsb ? fabsa*sqrt(1.0+std::pow(fabsb/fabsa,2)) :
	  (fabsb == 0.0 ? 0.0 : fabsb*sqrt(1.0+std::pow(fabsa/fabsb,2))));
}

