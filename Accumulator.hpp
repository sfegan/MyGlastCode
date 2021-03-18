//-*-mode:c++; mode:font-lock;-*-

#ifndef ACCUMULATOR_HPP
#define ACCUMULATOR_HPP

class Accumulator
{
  // Kahan summation algorithm - summation with correction
  // Reference: http://en.wikipedia.org/wiki/Kahan_summation_algorithm
public:
  Accumulator(): m_S(0), m_C(0) { }
  void reset(double s=0) { m_S = s; m_C = 0; }
  void add(double x)
  {
    /*volatile*/ double Y(x-m_C);
    /*volatile*/ double T(m_S+Y);
    m_C=(T-m_S)-Y;
    m_S=T;
  }
  double sum() const { return m_S; }
  double C() const { return m_C; }
private:
  double m_S;
  double m_C;
};

class TrapQuadIntegrator
{
public:
  TrapQuadIntegrator(double dx=1.0): 
    m_dx(dx), m_a(), m_have_one(false), m_x0(), m_xN() { }
  void reset(double dx=1.0) 
  { m_dx=dx; m_a.reset(); m_have_one=false; m_x0=0; m_xN=0; }
  void add(double x)
  {
    if(!m_have_one)m_x0=x;
    else { if(m_xN)m_a.add(m_xN); m_xN=x; }
  }
  double integral() const 
  { 
    Accumulator as(m_a); 
    as.add(0.5*m_x0); 
    as.add(0.5*m_xN);
    return as.sum()*m_dx;
  }
  double sum() const { return integral(); }
private:
  double m_dx;
  Accumulator m_a;
  bool m_have_one;
  double m_x0;
  double m_xN;
  
};

#endif
