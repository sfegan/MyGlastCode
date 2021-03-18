//-*-mode:c++; mode:font-lock;-*-

#ifndef FT2ROI_HPP
#define FT2ROI_HPP

#include <cmath>

#include "BLOBSerializer.hpp"
#include "FT2.hpp"
#include "GTI.hpp"
#include "Accumulator.hpp"
#include "PLIntegratedEA.hpp"

// ****************************************************************************
//
// FT2ROI: Integrate exposure over FT2 entries using effective area vs theta
//
// ****************************************************************************

class FT2ROI: public FITSVectorVisitor<FT2>
{
public:

  class Interval
  {
  public:
    Interval(): start(), stop() {}
    Interval(double t): start(t), stop(t) { }
    Interval(double _start, double _stop): start(_start), stop(_stop) { }
    double start;
    double stop;
    bool operator<(const Interval& o) const { return start<o.start; }
  };

  FT2ROI(const GTIRange& gti, double ra, double dec):
    FITSVectorVisitor<FT2>(), m_gti(gti), m_ft2gti(),
    m_ra(ra), m_dec(dec), m_zmax(M_PI), m_tmin(0), m_tmax(M_PI),
    m_rnlm(false) { }
  ~FT2ROI();

  void setZMaxCut(double zmax = M_PI) { m_zmax = zmax; }
  void setTMaxCut(double tmax = M_PI) { m_tmax = tmax; }
  void setTMinCut(double tmin = 0)    { m_tmin = tmin; }
  void setRequireNominalLATMode(bool rnlm = true) { m_rnlm = rnlm; }

  virtual void visitElement(unsigned irow, FT2& ft2);
  virtual void visitAcceptedInterval(unsigned irow, FT2& ft2);
  virtual void visitRejectedInterval(unsigned irow, FT2& ft2);

  bool isInFT2GTI(double t) const 
  {
    return m_ft2gti.findGTI(t)!=m_ft2gti.end(); 
  }

  const GTIRange& gti() const { return m_ft2gti; }

  bool serializeToBlob(BLOBSerializer& s) const;
  bool unserializeFromBlob(BLOBUnserializer& s);

protected:

  template<typename T> static 
  unsigned _findInterval(double t, const std::vector<T>& v)
  {
    typename std::vector<T>::const_iterator it = 
      std::upper_bound(v.begin(),v.end(),T(t));
    if(it==v.begin())it=v.end();
    else
      {
	it--;
	if((t<it->start)||(t>=it->stop))it=v.end();
      }
    return it-v.begin();
  }

  GTIRange        m_gti;
  GTIRange        m_ft2gti;
  double          m_ra;
  double          m_dec;
  double          m_zmax;
  double          m_tmin;
  double          m_tmax;
  bool            m_rnlm;
};

#endif // #ifndef FT2ROI_HPP
