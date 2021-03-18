//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAAlgebra.hpp
  Various linear algebra elements

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSAALGEBRA_HPP
#define VSAALGEBRA_HPP

#define VSA_ASSERT

#include<cmath>
#include<algorithm>
#include<vector>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<cassert>

namespace VERITAS
{
  namespace VSAAlgebra
  {

    // ========================================================================
    // 
    // Vec2D -- A two dimensional vector class
    //
    // ========================================================================

    class Vec2D
    {
    public:
      Vec2D(): fCx(), fCy() { /* nothing to see here */ }
      Vec2D(double cx, double cy): fCx(cx), fCy(cy) { /* nothing 2 C here */ }
      Vec2D(const Vec2D& v): fCx(v.fCx), fCy(v.fCy) { /* nothing 2 C here */ }

      // Getters
      double x() const { return fCx; }
      double y() const { return fCy; }
      double norm2() const { return fCx*fCx+fCy*fCy; }
      double norm() const { return std::sqrt(norm2()); }
      
      // Setters
      void set(double cx, double cy) { fCx=cx; fCy=cy; }
      void setNormalized(double cx, double cy) { set(cx,cy); normalize(); }
      inline void normalize();

      // Assignment Operators
      inline Vec2D& operator = (const Vec2D& v);  //! assignment
      inline Vec2D& operator += (const Vec2D& v); //! assignment addition
      inline Vec2D& operator -= (const Vec2D& v); //! assignment subtraction
      inline Vec2D& operator *= (double d); //! assignment scalar mult
      inline Vec2D& operator /= (double d); //! assignment scalar division

      // Non-assignment Operators
      inline Vec2D  operator + (const Vec2D& v) const; //! addition
      inline Vec2D  operator - (const Vec2D& v) const; //! subtraction
      inline double operator * (const Vec2D& v) const; //! scalar product
      inline Vec2D  operator - () const; //<negation

      // Miscellaneous
      void rotate(double theta);  //! rotates vector around axis
      void rotate(const Vec2D& point, double theta);
      void perp(Vec2D& p1) const;

      // Static Creation Functions
      static Vec2D make(double cx, double cy) { return Vec2D(cx,cy); }
      inline static Vec2D makeNormalized(double cx, double cy);
      inline static Vec2D makePolar(double theta, double r=1.0);
      
    private:
      double fCx;
      double fCy;
    };

    inline Vec2D operator * (double d, const Vec2D& v); //! left scalar mult.
    inline Vec2D operator * (const Vec2D& v, double d); //! right scalar mult
    inline Vec2D operator / (const Vec2D& v, double d); //! right scalar div

    // ========================================================================
    // 
    // Vec3D -- A three dimensional vector class
    //
    // ========================================================================

    enum Subspace2D { SS_12, SS_23, SS_31 };

    class RotationVec3D;

    class Vec3D
    {
    public:
      Vec3D(): fCx(), fCy(), fCz() { /* nothing to see here */ }
      Vec3D(double cx, double cy, double cz): fCx(cx), fCy(cy), fCz(cz) { }
      Vec3D(const Vec3D& v): fCx(v.fCx), fCy(v.fCy), fCz(v.fCz) { /* 02CH */ }

      // Getters
      double x() const { return fCx; }
      double y() const { return fCy; }
      double z() const { return fCz; }
      double a() const { return fCx; }
      double b() const { return fCy; }
      double c() const { return fCz; }
      double norm2() const { return fCx*fCx+fCy*fCy+fCz*fCz; }
      double norm() const { return std::sqrt(norm2()); }
      double theta() const { return std::atan2(std::sqrt(fCx*fCx+fCy*fCy),fCz);}
      double cosTheta() const { return fCz/norm(); }
      double phi() const { return std::atan2(fCy,fCx); }
      double lat() const { return std::atan2(fCz,std::sqrt(fCx*fCx+fCy*fCy)); }
      double lon() const { return phi(); }
      double lonPos() const { return fmod(phi()+2.0*M_PI,2.0*M_PI); }

      // Setters
      void set(double cx, double cy, double cz) { fCx=cx; fCy=cy; fCz=cz; }
      inline void setNormalized(double cx, double cy, double cz);
      inline void normalize();

      // Assignment Operators
      inline Vec3D& operator = (const Vec3D& v);  //! assignment
      inline Vec3D& operator += (const Vec3D& v); //! assignment addition
      inline Vec3D& operator -= (const Vec3D& v); //! assignment subtraction
      inline Vec3D& operator ^= (const Vec3D& v); //! assignment vector prod
      inline Vec3D& operator *= (double d); //! assignment scalar mult
      inline Vec3D& operator /= (double d); //! assignment scalar division
      
      Vec3D& operator &= (const Vec3D& r); //! assignment composition of rotations

      // Non-assignment Operators
      inline Vec3D  operator + (const Vec3D& v) const; //! addition
      inline Vec3D  operator - (const Vec3D& v) const; //! subtraction
      inline double operator * (const Vec3D& v) const; //! scalar product
      inline Vec3D  operator - () const; //! negation
      inline Vec3D  operator ^ (const Vec3D& v) const;  //! vector product
      inline Vec3D  operator & (const Vec3D& v) const;  //! composition of rotations

      // Miscellaneous
      void rotate(const Vec3D& axis);  //! rotates vector around axis
      void rotate(const RotationVec3D& axis);  //! rotates vector around axis
      void perp(Vec3D& p1, Vec3D& p2) const;
      void scatter(double sigma, double uniform_dev_1, double uniform_dev_2);

      // Static Creation Functions
      inline static Vec3D make(double cx, double cy, double cz);
      inline static Vec3D makeNormalized(double cx, double cy, double cz);
      inline static Vec3D makePolar(double theta, double phi, double r=1.0);
      inline static Vec3D makeLatLon(double lat, double lon, double r=1.0);

      static Vec3D makeRotation(const Vec3D& v1, const Vec3D& v2, 
				Subspace2D subspace = SS_12);
      inline static Vec3D makeRotationQuat(double w, 
					   double x, double y, double z);
      inline static Vec3D makeRotationQuat(double w, const Vec3D& v);

      inline static Vec3D makeRotationGalToJ2000()
      {
	//Vec3D rot(0,0,d2r(180-122.931918));
	//rot &= Vec3D(0,d2r(90.0-27.128251),0);
	//rot &= Vec3D(0,0,d2r(192.859481));
	return Vec3D(+1.174260663979527e+00,
		     -4.769204834242391e-01,
		     -1.699213212493570e+00);
      };
      
    protected:
      double fCx;
      double fCy;
      double fCz;
    };

    inline Vec3D operator * (double d, const Vec3D& v); //! left scalar mult.
    inline Vec3D operator * (const Vec3D& v, double d); //! right scalar mult
    inline Vec3D operator / (const Vec3D& v, double d); //! right scalar div

    // ========================================================================
    // 
    // RotationVec3D -- A three dimensional vector class for faster rotations
    //
    // ========================================================================

    class RotationVec3D: private Vec3D
    {
    public:
      inline RotationVec3D();
      inline RotationVec3D(double cx, double cy, double cz);
      inline RotationVec3D(double cx, double cy, double cz, 
			   double cos_theta, double sin_theta);
      inline RotationVec3D(const Vec3D& v);
      inline RotationVec3D(const RotationVec3D& rv);

      // Getters
      double x() const { return Vec3D::x(); }
      double y() const { return Vec3D::y(); }
      double z() const { return Vec3D::z(); }
      double a() const { return Vec3D::x(); }
      double b() const { return Vec3D::y(); }
      double c() const { return Vec3D::z(); }
      double cosTheta() const { return fCosTheta; }
      double sinTheta() const { return fSinTheta; }
      double theta() const { return std::atan2(fSinTheta,fCosTheta); }
      const Vec3D& axis() const { return *this; }

      //! Assignment composition of rotations
      RotationVec3D& operator &= (const RotationVec3D& r);
      inline RotationVec3D& operator &= (const Vec3D& r);

      //! Composition of rotations
      inline RotationVec3D  operator & (const RotationVec3D& v) const;
      inline RotationVec3D  operator & (const Vec3D& v) const;

      inline RotationVec3D  operator - () const; //! negation

      // Static Creation Functions
      static RotationVec3D makeRotation(const Vec3D& v1, const Vec3D& v2, 
					Subspace2D subspace = SS_12);
      inline static RotationVec3D makeRotationQuat(double w, 
						 double x, double y, double z);
      inline static RotationVec3D makeRotationQuat(double w, const Vec3D& v);

    private:
      double fCosTheta;
      double fSinTheta;
    };

    // ========================================================================
    // 
    // Ray3D -- Weighted 3D vector
    //
    // ========================================================================

    class Ray3D: public Vec3D
    {
    public:
      Ray3D(double weight, double cx, double cy, double cz)
	: Vec3D(cx,cy,cz), fWeight(weight) { }
      Ray3D(double weight, const Vec3D& vec)
	: Vec3D(vec), fWeight(weight) { }
      Ray3D(const Ray3D& ray)
	: Vec3D(ray), fWeight(ray.fWeight) { }
      
      double weight() const { return fWeight; }      
    private:
      double fWeight;
    };

    // ========================================================================
    // 
    // Symmetric2D -- A class for two dimensional real symmetric matrices
    //
    // ========================================================================

    struct Eigen2D
    {
      Eigen2D(): val(), vec() { }
      double val[2];
      Vec2D vec[2];
    };

    class Symmetric2D
    {
    public:
      Symmetric2D()
	: fA11(0), fA22(0), fA12(0) { /* nothing to see here */ }

      Symmetric2D(double a11, double a22, double a12) 
	: fA11(a11), fA22(a22), fA12(a12) { /* nothing to see here */ } 

      Symmetric2D(const Symmetric2D& m)
	: fA11(m.fA11), fA22(m.fA22), fA12(m.fA12) { }

      Symmetric2D(const Vec2D& v)
	: fA11(v.x()*v.x()), fA22(v.y()*v.y()), fA12(v.x()*v.y()) { }

      // Element access
      
      double a11() const { return fA11; }
      double a22() const { return fA22; }
      double a12() const { return fA12; }

      double xx() const { return fA11; }
      double yy() const { return fA22; }
      double xy() const { return fA12; }

      double aa() const { return fA11; }
      double bb() const { return fA22; }
      double ab() const { return fA12; }

      double maxAbsEl() const 
      {
	double max = fabs(fA11);
	double test;
	test = fabs(fA22); if(test>max)max=test;
	test = fabs(fA12); if(test>max)max=test;
	return max;
      }

      // Setters

      void set(double a11, double a22, double a12)  
      { fA11=a11; fA22=a22; fA12=a12; }

      void clear() { fA11=0; fA22=0; fA12=0; }

      // Invariants

      double trace() const { return fA11+fA22; }
      double det() const { return fA11*fA22 - fA12*fA12; }

      double invariant1() const { return trace(); }
      double invariant2() const { return det(); }

      // Inverse and eigenvalues

      inline Symmetric2D inverse() const;

      unsigned eigen(Eigen2D& e) const;
      unsigned eigenOld(Eigen2D& e) const;

      // Operators

      inline Symmetric2D& operator = (const Symmetric2D& m);
      inline Symmetric2D& operator += (const Symmetric2D& m);
      inline Symmetric2D& operator -= (const Symmetric2D& m);
      inline Symmetric2D& operator *= (double d);
      inline Symmetric2D& operator /= (double d);

      inline Symmetric2D& accumulate(double a11, double a22, double a12);

      inline Symmetric2D operator + (const Symmetric2D& m) const;
      inline Symmetric2D operator - (const Symmetric2D& m) const;

#ifndef SWIG
      const static Symmetric2D kronecker;
#endif

    private:
      double fA11;
      double fA22;
      double fA12;
    };

    inline Symmetric2D operator * (double d, const Symmetric2D& m);
    inline Symmetric2D operator * (const Symmetric2D& m, double d);
    inline Symmetric2D operator / (const Symmetric2D& m, double d);

    inline Vec2D operator * (const Symmetric2D& m, const Vec2D& v);
    inline Vec2D operator * (const Vec2D& v, const Symmetric2D& m);

    // ========================================================================
    // 
    // Symmetric3D -- A class for three dimensional real symmetric matrices
    //
    // ========================================================================
    
    struct Eigen3D
    {
      Eigen3D(): val(), vec() { }
      void clear() { val[0]=val[1]=val[2]=0; vec[0]=vec[1]=vec[2]=Vec3D(); }
      double val[3];
      Vec3D vec[3];
    };

    class Symmetric3D
    {
    public:
      Symmetric3D()
	: fA11(0), fA22(0), fA33(0), fA12(0), fA13(0), fA23(0) { }

      Symmetric3D(double a11, double a22, double a33, 
		  double a12, double a13, double a23)
	: fA11(a11), fA22(a22), fA33(a33), fA12(a12), fA13(a13), fA23(a23) { }

      Symmetric3D(const Symmetric3D& m)
	: fA11(m.fA11), fA22(m.fA22), fA33(m.fA33),
	  fA12(m.fA12), fA13(m.fA13), fA23(m.fA23) { }

      Symmetric3D(const Vec3D& v)
	: fA11(v.x()*v.x()), fA22(v.y()*v.y()), fA33(v.z()*v.z()), 
	  fA12(v.x()*v.y()), fA13(v.x()*v.z()), fA23(v.y()*v.z()) { }

      Symmetric3D(const Vec3D& v1, const Vec3D& v2)
	: fA11(2.0*v1.x()*v2.x()), 
	  fA22(2.0*v1.y()*v2.y()), 
	  fA33(2.0*v1.z()*v2.z()), 
	  fA12(v1.x()*v2.y()+v2.x()*v1.y()),
	  fA13(v1.x()*v2.z()+v2.x()*v1.z()),
	  fA23(v1.y()*v2.z()+v2.y()*v1.z()) { }
      

      // Element access
      
      double a11() const { return fA11; }
      double a22() const { return fA22; }
      double a33() const { return fA33; }
      double a12() const { return fA12; }
      double a21() const { return fA12; }
      double a13() const { return fA13; }
      double a31() const { return fA13; }
      double a23() const { return fA23; }
      double a32() const { return fA23; }

      double xx() const { return fA11; }
      double yy() const { return fA22; }
      double zz() const { return fA33; }
      double xy() const { return fA12; }
      double yx() const { return fA12; }
      double xz() const { return fA13; }
      double zx() const { return fA13; }
      double yz() const { return fA23; }
      double zy() const { return fA23; }

      double aa() const { return fA11; }
      double bb() const { return fA22; }
      double cc() const { return fA33; }
      double ab() const { return fA12; }
      double ba() const { return fA12; }
      double ac() const { return fA13; }
      double ca() const { return fA13; }
      double bc() const { return fA23; }
      double cb() const { return fA23; }

      double maxAbsDiagEl() const 
      {
	double max = fabs(fA11);
	double test;
	test = fabs(fA22); if(test>max)max=test;
	test = fabs(fA33); if(test>max)max=test;
	return max;
      }

      double maxAbsEl() const 
      {
	double max = fabs(fA11);
	double test;
	test = fabs(fA22); if(test>max)max=test;
	test = fabs(fA33); if(test>max)max=test;
	test = fabs(fA12); if(test>max)max=test;
	test = fabs(fA13); if(test>max)max=test;
	test = fabs(fA23); if(test>max)max=test;
	return max;
      }

      inline Symmetric2D subSpace(Subspace2D subspace) const;

      // Setters

      void set(double a11, double a22, double a33, 
	       double a12, double a13, double a23)
      { fA11=a11; fA22=a22; fA33=a33; fA12=a12; fA13=a12; fA23=a23; }

      void clear() { fA11=0; fA22=0; fA33=0; fA12=0; fA13=0; fA23=0; }

      // Invariants

      double trace() const { return fA11+fA22+fA33; }
      double det() const { return fA11*fA22*fA33 + 2.0*fA12*fA13*fA23 
	  - fA11*fA23*fA23 - fA22*fA13*fA13 - fA33*fA12*fA12; }

      double invariant1() const { return trace(); }
      double invariant2() const { return fA11*fA22 + fA11*fA33 + fA22*fA33
	  - fA12*fA12 - fA13*fA13 - fA23*fA23; }
      double invariant3() const { return det(); }

      // Inverse, eigenvalues and eigenvectors

      inline Symmetric3D inverse() const;

      unsigned eigen(Eigen3D& e) const { return /* eigenJacobi(e) */ eigenQL(e); }
      unsigned eigenQL(Eigen3D& e) const;
      unsigned eigenJacobi(Eigen3D& e) const;

      // Operators

      inline Symmetric3D& operator = (const Symmetric3D& m);
      inline Symmetric3D& operator += (const Symmetric3D& m);
      inline Symmetric3D& operator -= (const Symmetric3D& m);
      inline Symmetric3D& operator *= (double d);
      inline Symmetric3D& operator /= (double d);

      inline Symmetric3D& accumulate(double a11, double a22, double a33, 
				     double a12, double a13, double a23);

      inline Symmetric3D operator + (const Symmetric3D& m) const;
      inline Symmetric3D operator - (const Symmetric3D& m) const;

      // Others

      inline double traceProduct(const Symmetric3D& m) const;

#ifndef SWIG
      const static Symmetric3D kronecker;
#endif

    private:
      double fA11;
      double fA22;
      double fA33;
      double fA12;
      double fA13;
      double fA23;
    };

    inline Symmetric3D operator * (double d, const Symmetric3D& m);
    inline Symmetric3D operator * (const Symmetric3D& m, double d);
    inline Symmetric3D operator / (const Symmetric3D& m, double d);

    inline Vec3D operator * (const Symmetric3D& m, const Vec3D& v);
    inline Vec3D operator * (const Vec3D& v, const Symmetric3D& m);

    // ========================================================================
    // 
    // Orthogonal3D -- Three dimensional orthogonal transformations
    //
    // ========================================================================

    class Orthogonal3D
    {
    public:
      Orthogonal3D()
	: fRotation(), 
	  fAX(1), fAY(0), fAZ(0),
	  fBX(0), fBY(1), fBZ(0),
	  fCX(0), fCY(0), fCZ(1) { /* nothing to see here */ }
      
      Orthogonal3D(const Vec3D& rotation)
	: fRotation(rotation),
	  fAX(0), fAY(0), fAZ(0),
	  fBX(0), fBY(1), fBZ(0),
	  fCX(0), fCY(0), fCZ(0)
      { calcElements(); }
      
      Orthogonal3D(const Vec3D& v1, const Vec3D& v2, 
		   Subspace2D subspace = SS_12);
      
      // Element Access
      
      double ax() const { return fAX; }
      double ay() const { return fAY; }
      double az() const { return fAZ; }
      double bx() const { return fBX; }
      double by() const { return fBY; }
      double bz() const { return fBZ; }
      double cx() const { return fCX; }
      double cy() const { return fCY; }
      double cz() const { return fCZ; }

      const Vec3D& rotation() const { return fRotation; }

      // Orthogonal Transformations

      inline Orthogonal3D inverse() const;
      
      inline Vec3D transFwd(const Vec3D& v) const;
      inline Vec3D transBwd(const Vec3D& v) const;
      inline Symmetric3D transFwd(const Symmetric3D& s) const;
      inline Symmetric3D transBwd(const Symmetric3D& s) const;

      // Operators

      inline Orthogonal3D& operator *= (Orthogonal3D& o);

      inline Orthogonal3D operator - () const;
      inline Orthogonal3D operator * (const Orthogonal3D& o) const;

    private:
      void calcElements();

      Vec3D fRotation;
      double fAX;
      double fAY;
      double fAZ;
      double fBX;
      double fBY;
      double fBZ;
      double fCX;
      double fCY;
      double fCZ;
    };

    // ========================================================================
    // 
    // VecND: N-dimensional vector
    //
    // ========================================================================

    class VecND
    {
    public:
      VecND(): fNDim(), fCi() { /* nothing to see here */ }
      VecND(unsigned ndim, double c_all=double()): 
	fNDim(ndim), fCi(new double[fNDim]) { set(c_all); }
      VecND(unsigned ndim, const double* ci):
	fNDim(ndim), fCi(new double[fNDim]) { set(ci); }
      VecND(const std::vector<double>& ci):
	fNDim(ci.size()), fCi(new double[fNDim]) { set(&ci.front()); }
      VecND(const VecND& v): 
	fNDim(v.fNDim), fCi(new double[fNDim]) { set(v.fCi); }
      
      ~VecND() { delete[] fCi; }

      // Getters

      unsigned size() const { return fNDim; }
      unsigned ndim() const { return fNDim; }
      bool empty() const { return fNDim==0; }
      
      const double* data() const { return fCi; }

      double operator[] (unsigned idim) const { return fCi[idim]; }
      double operator() (unsigned idim) const { return fCi[idim]; }
      double at(unsigned idim) const { assert(idim<fNDim); return fCi[idim]; }

      double ci(unsigned idim) const { 
#ifdef VSA_ASSERT
	assert(idim<fNDim);
#endif // VSA_ASSERT	
	return fCi[idim]; 
      }
      
      // Setters

      VecND& operator=(const VecND& v) { set(v.fNDim, v.fCi); return *this; }

      void clear() { set(double()); }

      void resize(unsigned ndim) 
      { if(ndim>fNDim) 
	  { delete[] fCi; fCi = new double[ndim]; } fNDim = ndim; }

      double& operator[] (unsigned idim) { return fCi[idim]; }
      double& operator() (unsigned idim) { return fCi[idim]; }
      double& at(unsigned idim) { assert(idim<fNDim); return fCi[idim]; }

      void setEl(unsigned idim, const double ci) 
      {
#ifdef VSA_ASSERT
	assert(idim<fNDim);
#endif // VSA_ASSERT	
	fCi[idim] = ci;
      }

      void set(double c_all) { std::fill_n(fCi, fNDim, c_all); }
      void set(const double* ci) { std::copy(ci, ci+fNDim, fCi); }
      void set(unsigned ndim, const double* ci) { resize(ndim); set(ci); }
      void set(const std::vector<double>& ci) { set(ci.size(), &ci.front()); }

      void setNormalized(const double* ci) { set(ci); normalize(); }
      void setNormalized(unsigned ndim, const double* ci)
      { resize(ndim); setNormalized(ci); }
      void setNormalized(const std::vector<double>& ci)
      { setNormalized(ci.size(), &ci.front()); }

      inline void normalize();

      // Mathematical assignment Operators
      inline VecND& operator += (const VecND& v); //! assignment addition
      inline VecND& operator -= (const VecND& v); //! assignment subtraction
      inline VecND& operator *= (double d); //! assignment scalar mult
      inline VecND& operator /= (double d); //! assignment scalar division
      
      // Non-assignment Operators
      inline VecND  operator + (const VecND& v) const; //! addition
      inline VecND  operator - (const VecND& v) const; //! subtraction
      inline double operator * (const VecND& v) const; //! scalar product
      inline VecND  operator - () const; //! negation

    private:
      unsigned fNDim;
      double*  fCi;
    };

    inline VecND operator * (double d, const VecND& v); //! left scalar mult.
    inline VecND operator * (const VecND& v, double d); //! right scalar mult
    inline VecND operator / (const VecND& v, double d); //! right scalar div

    // ========================================================================
    // 
    // MatrixND: N-dimensional matrix
    //
    // ========================================================================

    class MatrixND
    {
    public:
      MatrixND(): fNRow(), fNCol(), fAij() { /* nothing to see here */ }

#if 0
      MatrixND(unsigned nrc, double a_all = double()): 
	fNRow(nrc), fNCol(nrc), fAij(new double[nel()]) { set(a_all); }
#endif
 
      MatrixND(unsigned nrow, unsigned ncol, double a_all = double()): 
	fNRow(nrow), fNCol(ncol), fAij(new double[nel()]) { set(a_all); }

#if 0
      MatrixND(unsigned nrc, const double* aij):
	fNRow(nrc), fNCol(nrc), fAij(new double[nel()]) { set(aij); }
#endif

      MatrixND(unsigned nrow, unsigned ncol, const double* aij):
	fNRow(nrow), fNCol(ncol), fAij(new double[nel()]) { set(aij); }

      MatrixND(const MatrixND& m):
	fNRow(m.fNRow), fNCol(m.fNCol), fAij(new double[nel()]){ set(m.fAij); }

      MatrixND(const VecND& v):
	fNRow(v.ndim()), fNCol(v.ndim()), fAij(new double[nel()])
      {
	for(unsigned irow=0;irow<fNRow;irow++)
	  for(unsigned icol=0;icol<fNCol;icol++)
	    *(fAij+ij(irow,icol)) = v[irow]*v[icol];
      }

      ~MatrixND() { delete[] fAij; }

      // Getters --------------------------------------------------------------

      unsigned nrow() const { return fNRow; }
      unsigned ncol() const { return fNCol; }
      unsigned ndim() const { return nrow(); }
      unsigned mdim() const { return ncol(); }
      bool empty() const { return fNRow==0 || fNCol==0; }

      bool isSquare() const { return fNRow==fNCol; }

      unsigned nel() const { return fNRow*fNCol; }
      unsigned ij(unsigned irow, unsigned icol) const
      { return irow*fNCol+icol; }
      
      const double* operator[] (unsigned irow) const { return &_el(irow,0); }

      const double operator() (unsigned irow, unsigned icol) const 
      { return _el(irow,icol); }

      const double at(unsigned irow, unsigned icol) const 
      { assert(irow<fNRow); assert(icol<fNCol); return _el(irow,icol); }
      
      const double aij(unsigned irow, unsigned icol) const 
      { 
#ifdef VSA_ASSERT
	assert(irow<fNRow);
	assert(icol<fNCol);
#endif // VSA_ASSERT	
	return _el(irow,icol); 
      }

      const double* row(unsigned irow) const { return &_el(irow,0); }

      // Setters --------------------------------------------------------------

      void clear() { set(double()); }

      void resize(unsigned nrow, unsigned ncol) 
      { if(nrow*ncol>nel()) { delete[] fAij; fAij = new double[nrow*ncol]; } 
	fNRow = nrow; fNCol = ncol; }

      MatrixND& operator = (const MatrixND& m) 
      { set(m.fNRow, m.fNCol, m.fAij); return *this; }

      double* operator[] (unsigned irow) { return &_el(irow,0); }
      double& operator() (unsigned irow, unsigned icol)
      { return _el(irow,icol); }
      double& at(unsigned irow, unsigned icol)
      { assert(irow<fNRow); assert(icol<fNCol); return _el(irow,icol); }
      
      void setEl(unsigned irow, unsigned icol, const double aij) 
      {
#ifdef VSA_ASSERT
	assert(irow<fNRow);
	assert(icol<fNCol);
#endif // VSA_ASSERT	
	_el(irow,icol) = aij;
      }
      
      void set(double a_all) { std::fill_n(fAij, nel(), a_all); }
      void set(const double* aij) 
      { std::copy(aij, aij+nel(), fAij); }
#if 0
      void set(unsigned nrc, const double* aij) 
      { resize(nrc,nrc); std::copy(aij, aij+nel(), fAij); }
#endif
      void set(unsigned nrow, unsigned ncol, const double* aij) 
      { resize(nrow,ncol); std::copy(aij, aij+nel(), fAij); }

      // Complex functions ----------------------------------------------------

      void invert()
      { assert(isSquare()); MatrixND b(fNRow,0,0.0); gaussJordan(*this,b); }
      void inverse(MatrixND &m) { m=*this; m.invert(); }
      MatrixND inverse() { MatrixND m; inverse(m); return m; }
      
      // Mathematical operators -----------------------------------------------

      inline MatrixND& operator += (const MatrixND& m);
      inline MatrixND& operator -= (const MatrixND& m);
      inline MatrixND& operator *= (const MatrixND& m);
      inline MatrixND& operator *= (double d);
      inline MatrixND& operator /= (double d);

      inline MatrixND operator + (const MatrixND& m) const;
      inline MatrixND operator - (const MatrixND& m) const;
      
      // Static functions -----------------------------------------------------

      static void gaussJordan(MatrixND& a, MatrixND& b);

    private:
      double& _el(unsigned irow, unsigned icol)
      { return *(fAij+ij(irow,icol)); }

      const double& _el(unsigned irow, unsigned icol) const
      { return *(fAij+ij(irow,icol)); }

      unsigned fNRow;
      unsigned fNCol;
      double* fAij;
    };
    
    inline MatrixND operator * (double d, const MatrixND& m);
    inline MatrixND operator * (const MatrixND& m, double d);
    inline MatrixND operator / (const MatrixND& m, double d);

    inline VecND operator * (const MatrixND& m, const VecND& v);
    inline VecND operator * (const VecND& v, const MatrixND& m);

    // ========================================================================
    // 
    // Vec2D: Definitions of longer inline functions
    //
    // ========================================================================
    
    inline void Vec2D::normalize()
    {
      double r = norm();
#ifdef VSA_ASSERT
      assert(r>0);
#endif // VSA_ASSERT	
      fCx/=r; 
      fCy/=r;
    }

    inline Vec2D& Vec2D::operator = (const Vec2D& v)
    {
      fCx=v.fCx;
      fCy=v.fCy;
      return *this;
    }

    inline Vec2D& Vec2D::operator += (const Vec2D& v)
    {
      fCx+=v.fCx;
      fCy+=v.fCy;
      return *this;
    }

    inline Vec2D& Vec2D::operator -= (const Vec2D& v)
    {
      fCx-=v.fCx;
      fCy-=v.fCy;
      return *this;
    }

    inline Vec2D& Vec2D::operator *= (double d)
    {
      fCx*=d;
      fCy*=d;
      return *this;
    }

    inline Vec2D& Vec2D::operator /= (double d)
    {
#ifdef VSA_ASSERT
      assert(fabs(d)>0);
#endif // VSA_ASSERT	
      fCx/=d;
      fCy/=d;
      return *this;
    }
    
    inline Vec2D Vec2D::operator + (const Vec2D& v) const
    {
      Vec2D o(*this);
      o+=v;
      return o;
    }

    inline Vec2D Vec2D::operator - (const Vec2D& v) const
    {
      Vec2D o(*this);
      o-=v;
      return o;
    }

    inline double Vec2D::operator * (const Vec2D& v) const
    {
      return fCx*v.fCx+fCy*v.fCy;
    }

    inline Vec2D Vec2D::operator - () const
    {
      return Vec2D(-fCx,-fCy);
    }
    
    inline Vec2D Vec2D::makeNormalized(double cx, double cy)
    { 
      Vec2D v(cx,cy);
      v.normalize();
      return v;
    }

    inline Vec2D Vec2D::makePolar(double theta, double r)
    {       
      return Vec2D(r*std::cos(theta), r*std::sin(theta));
    }

    inline Vec2D operator * (double d, const Vec2D& v)
    {
      Vec2D o(v);
      o*=d;
      return o;
    }

    inline Vec2D operator * (const Vec2D& v, double d)
    {
      return v*d;
    }

    inline Vec2D operator / (const Vec2D& v, double d)
    {
      Vec2D o(v);
      o/=d;
      return o;      
    }

    // ========================================================================
    //
    // Vec3D: Definitions of longer inline functions
    //
    // ========================================================================

    inline void Vec3D::setNormalized(double cx, double cy, double cz)
    {
      set(cx,cy,cz);
      normalize();
    }

    inline void Vec3D::normalize()
    {
      double r = norm();
#ifdef VSA_ASSERT
      assert(r>0);
#endif // VSA_ASSERT	
      fCx/=r; 
      fCy/=r;
      fCz/=r;
    }
    
    inline Vec3D& Vec3D::operator = (const Vec3D& v)
    {
      fCx=v.fCx;
      fCy=v.fCy;
      fCz=v.fCz;
      return *this;
    }

    inline Vec3D& Vec3D::operator += (const Vec3D& v)
    {
      fCx+=v.fCx;
      fCy+=v.fCy;
      fCz+=v.fCz;
      return *this;      
    }

    inline Vec3D& Vec3D::operator -= (const Vec3D& v)
    {
      fCx-=v.fCx;
      fCy-=v.fCy;
      fCz-=v.fCz;
      return *this;      
    }

    inline Vec3D& Vec3D::operator ^= (const Vec3D& v)
    {
      set(fCy*v.fCz - fCz*v.fCy, fCz*v.fCx - fCx*v.fCz, fCx*v.fCy - fCy*v.fCx);
      return *this;
    }

    inline Vec3D& Vec3D::operator *= (double d)
    {
      fCx*=d;
      fCy*=d;
      fCz*=d;
      return *this;
    }

    inline Vec3D& Vec3D::operator /= (double d)
    {
#ifdef VSA_ASSERT
      assert(fabs(d)>0);
#endif // VSA_ASSERT	
      fCx/=d;
      fCy/=d;
      fCz/=d;
      return *this;
    }
    
    // Non-assignment Operators
    inline Vec3D Vec3D::operator + (const Vec3D& v) const
    {
      Vec3D o(*this);
      o+=v;
      return o;
    }

    inline Vec3D Vec3D::operator - (const Vec3D& v) const
    {
      Vec3D o(*this);
      o-=v;
      return o;
    }

    inline double Vec3D::operator * (const Vec3D& v) const
    {
      return fCx*v.fCx+fCy*v.fCy+fCz*v.fCz;      
    }

    inline Vec3D Vec3D::operator - () const
    {
      return Vec3D(-fCx,-fCy,-fCz);      
    }

    inline Vec3D Vec3D::operator ^ (const Vec3D& v) const
    {
      Vec3D o(*this);
      o ^= v;
      return o;
    }

    inline Vec3D Vec3D::operator & (const Vec3D& v) const
    {
      Vec3D o(*this);
      o &= v;
      return o;
    }
    
    inline Vec3D Vec3D::make(double cx, double cy, double cz) 
    {
      return Vec3D(cx,cy,cz);
    }    
    
    inline Vec3D Vec3D::makeNormalized(double cx, double cy, double cz)
    {
      Vec3D v(cx,cy,cz);
      v.normalize();
      return v;
    }

    inline Vec3D Vec3D::makePolar(double theta, double phi, double r)
    { 
      const double rs = r*std::sin(theta);
      return Vec3D(rs*std::cos(phi), rs*std::sin(phi), r*std::cos(theta));
    }

    inline Vec3D Vec3D::makeLatLon(double lat, double lon, double r)
    {
      return makePolar(M_PI_2-lat,lon,r);
    }

    inline Vec3D 
    Vec3D::makeRotationQuat(double w, double x, double y, double z)
    {
      Vec3D v(x,y,z);
      const double st2 = v.norm();
      v *= 2.0*std::atan2(st2,w)/st2;
      return v;
    }
     
    inline Vec3D 
    Vec3D::makeRotationQuat(double w, const Vec3D& v)
    {
      return makeRotationQuat(w, v.x(), v.y(), v.z());
    }

    inline Vec3D operator * (double d, const Vec3D& v)
    {
      Vec3D o(v);
      o*=d;
      return o;
    }

    inline Vec3D operator * (const Vec3D& v, double d)
    {
      return d*v;
    }

    inline Vec3D operator / (const Vec3D& v, double d)
    {
      Vec3D o(v);
      o/=d;
      return o;
    }

    // ========================================================================
    //
    // RotationRay3D: Definitions of longer inline functions
    //
    // ========================================================================

    inline RotationVec3D::RotationVec3D(): 
      Vec3D(), fCosTheta(), fSinTheta()
    {
      // nothing to see here    
    }

    inline RotationVec3D::RotationVec3D(double cx, double cy, double cz): 
      Vec3D(cx, cy, cz), fCosTheta(), fSinTheta() 
    { 
      double norm = std::sqrt(fCx*fCx+fCy*fCy+fCz*fCz);
      double norm_inv = 1.0/norm;
      fCx *= norm_inv;
      fCy *= norm_inv;
      fCz *= norm_inv;
      fCosTheta = std::cos(norm);
      fSinTheta = std::sin(norm);
    }

    inline RotationVec3D::RotationVec3D(double cx, double cy, double cz, 
					double cos_theta, double sin_theta):
      Vec3D(cx, cy, cz), fCosTheta(cos_theta), fSinTheta(sin_theta)
    {
      // nothing to see here
    }

    inline RotationVec3D::RotationVec3D(const Vec3D& v):
      Vec3D(v), fCosTheta(), fSinTheta() 
    { 
      double norm = std::sqrt(fCx*fCx+fCy*fCy+fCz*fCz);
      double norm_inv = 1.0/norm;
      fCx *= norm_inv;
      fCy *= norm_inv;
      fCz *= norm_inv;
      fCosTheta = std::cos(norm);
      fSinTheta = std::sin(norm);
    }

    inline RotationVec3D::RotationVec3D(const RotationVec3D& rv):
      Vec3D(rv), fCosTheta(rv.fCosTheta), fSinTheta(rv.fSinTheta)
    {
      // nothing to see here    
    }

    inline RotationVec3D& RotationVec3D::operator &= (const Vec3D& r)
    {
      Vec3D t(axis()*theta());
      t &= r;
      *this = t;
      return *this;
    }

    inline RotationVec3D RotationVec3D::operator & (const RotationVec3D& v) const
    {
      RotationVec3D o(*this);
      o &= v;
      return o;
    }

    inline RotationVec3D RotationVec3D::operator & (const Vec3D& v) const
    {
      RotationVec3D o(*this);
      o &= v;
      return o;
    }

    inline RotationVec3D RotationVec3D::operator - () const
    {
      return RotationVec3D(fCx,fCy,fCz,fCosTheta,-fSinTheta);
    }

    inline RotationVec3D 
    RotationVec3D::makeRotationQuat(double w, double x, double y, double z)
    {
      const double st2 = std::sqrt(x*x+y*y+z*z);
      const double ct2 = w;
      const double norm_inv = 1.0/st2;
      x *= norm_inv;
      y *= norm_inv;
      z *= norm_inv;
      return RotationVec3D(x,y,z,ct2*ct2-st2*st2,2.0*st2*ct2);
    }
     
    inline RotationVec3D 
    RotationVec3D::makeRotationQuat(double w, const Vec3D& v)
    {
      return makeRotationQuat(w, v.x(), v.y(), v.z());
    }

    // ========================================================================
    //
    // Ray3D: Definitions of longer inline functions
    //
    // ========================================================================


    // ========================================================================
    //
    // Symmetric2D: Definitions of longer inline functions
    //
    // ========================================================================
    
    inline Symmetric2D Symmetric2D::inverse() const
    {
      double d = det();
#ifdef VSA_ASSERT
      assert(fabs(d)>0);
#endif // VSA_ASSERT
      return Symmetric2D(fA22, fA11, -fA12)/d;
    }

    inline Symmetric2D& Symmetric2D::operator = (const Symmetric2D& m)
    {
      fA11 = m.fA11;
      fA22 = m.fA22;
      fA12 = m.fA12;
      return *this;
    }

    inline Symmetric2D& Symmetric2D::operator += (const Symmetric2D& m)
    {
      fA11 += m.fA11;
      fA22 += m.fA22;
      fA12 += m.fA12;
      return *this;
    }

    inline Symmetric2D& Symmetric2D::operator -= (const Symmetric2D& m)
    {
      fA11 -= m.fA11;
      fA22 -= m.fA22;
      fA12 -= m.fA12;
      return *this;
    }

    inline Symmetric2D& Symmetric2D::operator *= (double d)
    {
      fA11 *= d;
      fA22 *= d;
      fA12 *= d;
      return *this;
    }

    inline Symmetric2D& Symmetric2D::operator /= (double d)
    {
      fA11 /= d;
      fA22 /= d;
      fA12 /= d;
      return *this;
    }
    
    inline Symmetric2D& 
    Symmetric2D::accumulate(double a11, double a22, double a12)
    {
      fA11 += a11;
      fA22 += a22;
      fA12 += a12;
      return *this;
    }
    
    inline Symmetric2D Symmetric2D::operator + (const Symmetric2D& m) const
    {
      Symmetric2D r(*this); r+=m;
      return r;
    }

    inline Symmetric2D Symmetric2D::operator - (const Symmetric2D& m) const
    {
      Symmetric2D r(*this); r-=m;
      return r;
    }
    
    inline Symmetric2D operator * (double d, const Symmetric2D& m)
    {
      Symmetric2D r(m); r*=d;
      return r;
    }

    inline Symmetric2D operator * (const Symmetric2D& m, double d)
    {
      Symmetric2D r(m); r*=d;
      return r;
    }

    inline Symmetric2D operator / (const Symmetric2D& m, double d)
    {
      Symmetric2D r(m); r/=d;
      return r;
    }

    inline Vec2D operator * (const Symmetric2D& m, const Vec2D& v)
    {
      return Vec2D(m.a11()*v.x() + m.a12()*v.y(),
		   m.a12()*v.x() + m.a22()*v.y());
    }

    inline Vec2D operator * (const Vec2D& v, const Symmetric2D& m)
    {
      return Vec2D(m.a11()*v.x() + m.a12()*v.y(),
		   m.a12()*v.x() + m.a22()*v.y());
    }

    // ========================================================================
    //
    // Symmetric3D: Definitions of longer inline functions
    //
    // ========================================================================

    inline Symmetric2D Symmetric3D::subSpace(Subspace2D subspace) const
    {
      switch(subspace)
	{
	case SS_12: return Symmetric2D(fA11,fA22,fA12);
	case SS_23: return Symmetric2D(fA22,fA33,fA23);
	case SS_31: return Symmetric2D(fA33,fA11,fA13);
	}
      assert(0);
    }

    inline Symmetric3D Symmetric3D::inverse() const
    {
      double d = det();
#ifdef VSA_ASSERT
      assert(fabs(d)>0);
#endif // VSA_ASSERT
      return Symmetric3D(fA22*fA33-fA23*fA23, fA11*fA33-fA13*fA13,
			 fA11*fA22-fA12*fA23, fA23*fA13-fA12*fA33,
			 fA12*fA23-fA22*fA13, fA13*fA12-fA11*fA23)/d;
    }

    inline Symmetric3D& Symmetric3D::operator = (const Symmetric3D& m)
    {
      fA11 = m.fA11; fA22 = m.fA22; fA33 = m.fA33;
      fA12 = m.fA12; fA13 = m.fA13; fA23 = m.fA23;
      return *this;
    }

    inline Symmetric3D& Symmetric3D::operator += (const Symmetric3D& m)
    {
      fA11 += m.fA11; fA22 += m.fA22; fA33 += m.fA33;
      fA12 += m.fA12; fA13 += m.fA13; fA23 += m.fA23;
      return *this;
    }

    inline Symmetric3D& Symmetric3D::operator -= (const Symmetric3D& m)
    {
      fA11 -= m.fA11; fA22 -= m.fA22; fA33 -= m.fA33;
      fA12 -= m.fA12; fA13 -= m.fA13; fA23 -= m.fA23;
      return *this;
    }
    
    inline Symmetric3D& Symmetric3D::operator *= (double d)
    {
      fA11 *= d; fA22 *= d; fA33 *= d;
      fA12 *= d; fA13 *= d; fA23 *= d;
      return *this;
    }

    inline Symmetric3D& Symmetric3D::operator /= (double d)
    {
      fA11 /= d; fA22 /= d; fA33 /= d;
      fA12 /= d; fA13 /= d; fA23 /= d;
      return *this;
    }

    inline Symmetric3D& 
    Symmetric3D::accumulate(double a11, double a22, double a33, 
			    double a12, double a13, double a23)
    {
      fA11 += a11; fA22 += a22; fA33 += a33;
      fA12 += a12; fA13 += a13; fA23 += a23;
      return *this;
    }
    
    inline Symmetric3D Symmetric3D::operator + (const Symmetric3D& m) const
    {
      Symmetric3D r(*this); r+=m;
      return r;
    }

    inline double Symmetric3D::traceProduct(const Symmetric3D& m) const
    {
      return fA11*m.fA11 + fA22*m.fA22 + fA33*m.fA33
	+ 2.0 * ( fA12*m.fA12 + fA13*m.fA13 + fA23*m.fA23 );
    }

    inline Symmetric3D Symmetric3D::operator - (const Symmetric3D& m) const
    {
      Symmetric3D r(*this); r-=m;
      return r;
    }

    inline Symmetric3D operator * (double d, const Symmetric3D& m)
    {
      Symmetric3D r(m); r*=d;
      return r;
    }

    inline Symmetric3D operator * (const Symmetric3D& m, double d)
    {
      Symmetric3D r(m); r*=d;
      return r;      
    }

    inline Symmetric3D operator / (const Symmetric3D& m, double d)
    {
      Symmetric3D r(m); r/=d;
      return r;
    }

    inline Vec3D operator * (const Symmetric3D& m, const Vec3D& v)
    {
      return Vec3D(m.a11()*v.x() + m.a12()*v.y() + m.a13()*v.z(),
		   m.a12()*v.x() + m.a22()*v.y() + m.a23()*v.z(),
		   m.a13()*v.x() + m.a23()*v.y() + m.a33()*v.z());
    }

    inline Vec3D operator * (const Vec3D& v, const Symmetric3D& m)
    {
      return Vec3D(m.a11()*v.x() + m.a12()*v.y() + m.a13()*v.z(),
		   m.a12()*v.x() + m.a22()*v.y() + m.a23()*v.z(),
		   m.a13()*v.x() + m.a23()*v.y() + m.a33()*v.z());
    }

    // ========================================================================
    //
    // Orthogonal3D: Definitions of longer inline functions
    //
    // ========================================================================

    inline Orthogonal3D Orthogonal3D::inverse() const
    {
      return Orthogonal3D(-fRotation);
    }

    inline Vec3D Orthogonal3D::transFwd(const Vec3D& v) const
    {
      const double va = fAX*v.x() + fAY*v.y() + fAZ*v.z();
      const double vb = fBX*v.x() + fBY*v.y() + fBZ*v.z();
      const double vc = fCX*v.x() + fCY*v.y() + fCZ*v.z();
      return Vec3D(va,vb,vc);
    }

    inline Vec3D Orthogonal3D::transBwd(const Vec3D& v) const
    {
      const double vx = v.a()*fAX + v.b()*fBX + v.c()*fCX;
      const double vy = v.a()*fAY + v.b()*fBY + v.c()*fCY;
      const double vz = v.a()*fAZ + v.b()*fBZ + v.c()*fCZ;
      return Vec3D(vx,vy,vz);
    }

    inline Symmetric3D Orthogonal3D::transFwd(const Symmetric3D& s) const
    {
      const double sax = fAX*s.xx() + fAY*s.yx() + fAZ*s.zx();
      const double say = fAX*s.xy() + fAY*s.yy() + fAZ*s.zy();
      const double saz = fAX*s.xz() + fAY*s.yz() + fAZ*s.zz();
      const double sbx = fBX*s.xx() + fBY*s.yx() + fBZ*s.zx();
      const double sby = fBX*s.xy() + fBY*s.yy() + fBZ*s.zy();
      const double sbz = fBX*s.xz() + fBY*s.yz() + fBZ*s.zz();
      const double scx = fCX*s.xx() + fCY*s.yx() + fCZ*s.zx();
      const double scy = fCX*s.xy() + fCY*s.yy() + fCZ*s.zy();
      const double scz = fCX*s.xz() + fCY*s.yz() + fCZ*s.zz();

      const double saa = sax*fAX + say*fAY + saz*fAZ;
      const double sab = sax*fBX + say*fBY + saz*fBZ;
      const double sac = sax*fCX + say*fCY + saz*fCZ;
      //const double sba = sbx*fAX + sby*fAY + sbz*fAZ;
      const double sbb = sbx*fBX + sby*fBY + sbz*fBZ;
      const double sbc = sbx*fCX + sby*fCY + sbz*fCZ;
      //const double sca = scx*fAX + scy*fAY + scz*fAZ;
      //const double scb = scx*fBX + scy*fBY + scz*fBZ;
      const double scc = scx*fCX + scy*fCY + scz*fCZ;

      return Symmetric3D(saa,sbb,scc,sab,sac,sbc);
    }

    inline Symmetric3D Orthogonal3D::transBwd(const Symmetric3D& s) const
    {
      const double sxa = fAX*s.aa() + fBX*s.ba() + fCX*s.ca();
      const double sxb = fAX*s.ab() + fBX*s.bb() + fCX*s.cb();
      const double sxc = fAX*s.ac() + fBX*s.bc() + fCX*s.cc();
      const double sya = fAY*s.aa() + fBY*s.ba() + fCY*s.ca();
      const double syb = fAY*s.ab() + fBY*s.bb() + fCY*s.cb();
      const double syc = fAY*s.ac() + fBY*s.bc() + fCY*s.cc();
      const double sza = fAZ*s.aa() + fBZ*s.ba() + fCZ*s.ca();
      const double szb = fAZ*s.ab() + fBZ*s.bb() + fCZ*s.cb();
      const double szc = fAZ*s.ac() + fBZ*s.bc() + fCZ*s.cc();

      const double sxx = sxa*fAX + sxb*fBX + sxc*fCX;
      const double sxy = sxa*fAY + sxb*fBY + sxc*fCY;
      const double sxz = sxa*fAZ + sxb*fBZ + sxc*fCZ;
      //const double syx = sya*fAX + syb*fBX + syc*fCX;
      const double syy = sya*fAY + syb*fBY + syc*fCY;
      const double syz = sya*fAZ + syb*fBZ + syc*fCZ;
      //const double szx = sza*fAX + szb*fBX + szc*fCX;
      //const double szy = sza*fAY + szb*fBY + szc*fCY;
      const double szz = sza*fAZ + szb*fBZ + szc*fCZ;

      return Symmetric3D(sxx,syy,szz,sxy,sxz,syz);
    }
    
    inline Orthogonal3D& Orthogonal3D::operator *= (Orthogonal3D& o)
    {
      // Rotation composition is opposite in sense to matrix multiplication
      fRotation = o.fRotation&fRotation;
      calcElements();
      return *this;
    }
    
    inline Orthogonal3D Orthogonal3D::operator - () const
    {
      return Orthogonal3D(-fRotation);
    }

    inline Orthogonal3D Orthogonal3D::operator * (const Orthogonal3D& o) const
    {
      return Orthogonal3D(o.fRotation&fRotation);
    }
   
    // ========================================================================
    //
    // VecND: Definitions of longer inline functions
    //
    // ========================================================================

    inline void VecND::normalize()
    {
      double x = 0;
      for(unsigned idim=0;idim<fNDim;idim++)x += fCi[idim]*fCi[idim];
      x = 1.0/std::sqrt(x);
      for(unsigned idim=0;idim<fNDim;idim++)fCi[idim] *= x;
    }

    inline VecND& VecND::operator += (const VecND& v)
    {
#ifdef VSA_ASSERT
      assert(fNDim==v.fNDim);
#endif // VSA_ASSERT	
      for(unsigned idim=0;idim<fNDim;idim++)fCi[idim] += v.fCi[idim];
      return *this;
    }

    inline VecND& VecND::operator -= (const VecND& v)
    {
#ifdef VSA_ASSERT
      assert(fNDim==v.fNDim);
#endif // VSA_ASSERT	
      for(unsigned idim=0;idim<fNDim;idim++)fCi[idim] -= v.fCi[idim];
      return *this;
    }

    inline VecND& VecND::operator *= (double d)
    {
      for(unsigned idim=0;idim<fNDim;idim++)fCi[idim] *= d;
      return *this;
    }

    inline VecND& VecND::operator /= (double d)
    {
      d = 1.0/d;
      for(unsigned idim=0;idim<fNDim;idim++)fCi[idim] *= d;
      return *this;
    }
    
    inline VecND VecND::operator + (const VecND& v) const
    {
      VecND vret(*this); return vret += v;
    }

    inline VecND VecND::operator - (const VecND& v) const
    {
      VecND vret(*this); return vret -= v;
    }

    inline double VecND::operator * (const VecND& v) const
    {
#ifdef VSA_ASSERT
      assert(fNDim==v.fNDim);
#endif // VSA_ASSERT	
      double x = 0;
      for(unsigned idim=0;idim<fNDim;idim++)x = fCi[idim]*v.fCi[idim];
      return x;
    }

    inline VecND VecND::operator - () const
    {
      VecND vret(fNDim);
      for(unsigned idim=0;idim<fNDim;idim++)vret.fCi[idim] = -fCi[idim];
      return vret;
    }

    inline VecND operator * (double d, const VecND& v)
    {
      VecND vret(v); return vret *= d;
    }
    
    inline VecND operator * (const VecND& v, double d)
    {
      VecND vret(v); return vret *= d;
    }

    inline VecND operator / (const VecND& v, double d)
    {
      VecND vret(v); return vret /= d;
    }

    // ========================================================================
    //
    // MatrixND: Definitions of longer inline functions
    //
    // ========================================================================

    inline MatrixND& MatrixND::operator += (const MatrixND& m)
    {
#ifdef VSA_ASSERT
      assert(fNRow == m.fNRow);
      assert(fNCol == m.fNCol);
#endif // VSA_ASSERT	
      unsigned _nel = nel();
      for(unsigned iel=0;iel<_nel;iel++)fAij[iel] += m.fAij[iel];
      return *this;
    }

    inline MatrixND& MatrixND::operator -= (const MatrixND& m)
    {
#ifdef VSA_ASSERT
      assert(fNRow == m.fNRow);
      assert(fNCol == m.fNCol);
#endif // VSA_ASSERT	
      unsigned _nel = nel();
      for(unsigned iel=0;iel<_nel;iel++)fAij[iel] -= m.fAij[iel];
      return *this;      
    }

    inline MatrixND& MatrixND::operator *= (const MatrixND& m)
    {
#ifdef VSA_ASSERT
      assert(isSquare());
      assert(m.isSquare());
      assert(fNCol == m.fNRow);
#endif // VSA_ASSERT	
      VecND col(fNCol);
      for(unsigned irow=0;irow<fNRow;irow++)
	{
	  col.set(row(irow));
	  for(unsigned icol=0;icol<fNCol;icol++)
	    {
	      double x = 0;
	      for(unsigned jrow=0;jrow<m.fNRow;jrow++)
		x += col[jrow]*m(jrow,icol);
	      _el(irow,icol) = x;
	    }
	}
      return *this;      
    }

    inline MatrixND& MatrixND::operator *= (double d)
    {
      unsigned _nel = nel();
      for(unsigned iel=0;iel<_nel;iel++)fAij[iel] *= d;
      return *this;      
    }

    inline MatrixND& MatrixND::operator /= (double d)
    {
      d = 1.0/d;
      unsigned _nel = nel();
      for(unsigned iel=0;iel<_nel;iel++)fAij[iel] *= d;
      return *this;      
    }
    
    inline MatrixND MatrixND::operator + (const MatrixND& m) const
    {
      MatrixND mret(*this); return mret += m;
    }

    inline MatrixND MatrixND::operator - (const MatrixND& m) const
    {
      MatrixND mret(*this); return mret -= m;
    }

    inline MatrixND operator * (double d, const MatrixND& m)
    {
      MatrixND mret(m); return mret *= d;
    }

    inline MatrixND operator * (const MatrixND& m, double d)
    {
      MatrixND mret(m); return mret *= d;
    }

    inline MatrixND operator / (const MatrixND& m, double d)
    {
      MatrixND mret(m); return mret /= d;
    }

    inline VecND operator * (const MatrixND& m, const VecND& v)
    {
      unsigned ncol = m.ncol();
      unsigned nrow = m.nrow();
#ifdef VSA_ASSERT
      assert(v.ndim()==ncol);
#endif // VSA_ASSERT	
      VecND vret(nrow);
      for(unsigned irow=0;irow<nrow;irow++)
	{
	  const double* row = m[irow];
	  double x = 0;
	  for(unsigned icol=0;icol<ncol;icol++)
	    x += row[icol]*v[icol];
	  vret[irow] = x;
	}
      return vret;
    }

    inline VecND operator * (const VecND& v, const MatrixND& m)
    {
      unsigned ncol = m.ncol();
      unsigned nrow = m.nrow();
#ifdef VSA_ASSERT
      assert(v.ndim()==nrow);
#endif // VSA_ASSERT	
      VecND vret(ncol);
      for(unsigned irow=0;irow<nrow;irow++)
	{
	  const double* row = m[irow];
	  double x = v[irow];
	  for(unsigned icol=0;icol<ncol;icol++)
	    vret[icol] += x*row[icol];
	}
      return vret;
    }

  }

}

inline std::ostream& operator<< (std::ostream& stream, 
				 const VERITAS::VSAAlgebra::Vec2D& v)
{
  std::ostringstream _stream;
  _stream << std::showpos << std::scientific << std::setprecision(8) 
	  << v.x() << ' ' << v.y();
  stream << _stream.str();
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, 
				 const VERITAS::VSAAlgebra::Vec3D& v)
{
  std::ostringstream _stream;
  _stream << std::showpos << std::scientific << std::setprecision(8) 
	  << v.x() << ' ' << v.y() << ' ' << v.z();
  stream << _stream.str();
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, 
				 const VERITAS::VSAAlgebra::Ray3D& v)
{
  std::ostringstream _stream;
  _stream << std::showpos << std::scientific << std::setprecision(8) 
	  << v.weight() << ' ' 
	  << static_cast<const VERITAS::VSAAlgebra::Vec3D&>(v);
  stream << _stream.str();
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, 
				 const VERITAS::VSAAlgebra::Symmetric2D& v)
{
  std::ostringstream _stream;
  _stream << std::showpos << std::scientific << std::setprecision(8) 
	  << v.xx() << ' ' << v.xy() << std::endl
	  << v.xy() << ' ' << v.yy();
  stream << _stream.str();
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, 
				 const VERITAS::VSAAlgebra::Symmetric3D& v)
{
  std::ostringstream _stream;
  _stream << std::showpos << std::scientific << std::setprecision(8) 
	  << v.xx() << ' ' << v.xy() << ' ' << v.xz() << std::endl
	  << v.xy() << ' ' << v.yy() << ' ' << v.yz() << std::endl
	  << v.xz() << ' ' << v.yz() << ' ' << v.zz();
  stream << _stream.str();
  return stream;
}

inline std::ostream& operator<< (std::ostream& stream, 
				 const VERITAS::VSAAlgebra::Orthogonal3D& v)
{
  std::ostringstream _stream;
  _stream << std::showpos << std::scientific << std::setprecision(8)
	  << v.ax() << ' ' << v.ay() << ' ' << v.az() << std::endl
	  << v.bx() << ' ' << v.by() << ' ' << v.bz() << std::endl
	  << v.cx() << ' ' << v.cy() << ' ' << v.cz();
  stream << _stream.str();
  return stream;
}

#endif // VSAALGEBRA_HP
