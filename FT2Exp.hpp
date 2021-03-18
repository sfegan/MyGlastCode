//-*-mode:c++; mode:font-lock;-*-

#ifndef FT2EXP_HPP
#define FT2EXP_HPP

#include "FT2ROI.hpp"
#include "Accumulator.hpp"
#include "PLIntegratedEA.hpp"

// ****************************************************************************
//
// FT2Exp: Integrate exposure over FT2 entries using effective area vs theta
//
// ****************************************************************************

class FT2Exp: public FT2ROI
{
public:

  class TimeInterval
  {
  public:
    TimeInterval(): 
      start(), stop(), dlivetime(), livetime_int_stop(), costheta(), phi() {}
    TimeInterval(double t): 
      start(t), stop(t), dlivetime(), livetime_int_stop(), costheta(), phi() { }
    TimeInterval(double _start, double _stop, double _dlivetime, 
		 double _livetime_int_stop, double _costheta, double _phi): 
      start(_start), stop(_stop), dlivetime(_dlivetime), 
      livetime_int_stop(_livetime_int_stop), costheta(_costheta), phi(_phi) { }
    union {
      double start;
      double clock_time_start;
    };
    union {
      double stop;
      double clock_time_stop;
    };
    double dlivetime;
    double livetime_int_stop;
    double costheta;
    double phi;
    bool operator<(const TimeInterval& o) const { return start<o.start; }

    bool serializeToBlob(BLOBSerializer& s) const
    { return s.serialize(clock_time_start) && s.serialize(clock_time_stop)
	&& s.serialize(dlivetime) && s.serialize(livetime_int_stop)
	&& s.serialize(costheta) && s.serialize(phi); }
    bool unserializeFromBlob(BLOBUnserializer& s)
    { return s.unserialize(clock_time_start) && s.unserialize(clock_time_stop)
	&& s.unserialize(dlivetime) && s.unserialize(livetime_int_stop)
	&& s.unserialize(costheta) && s.unserialize(phi); }
  };
  
  class ExposureInterval
  {
  public:
    ExposureInterval(): start(), stop(), aeff(), dexp() {}
    ExposureInterval(double t): start(t), stop(t), aeff(), dexp() { }
    ExposureInterval(double _start, double _stop, double _aeff, double _dexp): 
      start(_start), stop(_stop), aeff(_aeff), dexp(_dexp) { }
    union {
      double start;
      double exp_int_start;
    };
    union {
      double stop;
      double exp_int_stop;
    };
    double aeff;
    double dexp;
    bool operator<(const ExposureInterval& o) const { return start<o.start; }
    bool serializeToBlob(BLOBSerializer& s) const
    { return s.serialize(exp_int_start) && s.serialize(exp_int_stop)
	&& s.serialize(aeff) && s.serialize(dexp); }
    bool unserializeFromBlob(BLOBUnserializer& s)
    { return s.unserialize(exp_int_start) && s.unserialize(exp_int_stop)
	&& s.unserialize(aeff) && s.unserialize(dexp); }
  };

  FT2Exp(const std::vector<PLIntegratedEA*>& ea, 
	 const GTIRange& gti, double ra, double dec, double lc_period_day=1.0): 
    FT2ROI(gti,ra,dec), m_t0set(false), m_t0(), 
    m_tint(), m_eint(ea.size()), m_ea(ea), m_tacc(), m_eacc(ea.size()), 
    m_merged_exposures(false), m_lcperiod(lc_period_day*86400.0), m_lcbin(), 
    m_lcacc(), m_lcexp() 
  {
    setRequireNominalLATMode();
  }
  ~FT2Exp();

  //virtual void visitElement(unsigned irow, FT2& ft2);
  virtual void visitAcceptedInterval(unsigned irow, FT2& ft2);
  virtual void leaveFileVector(); 

  void mergeExposures();
  bool hasMergedExposures() const { return m_merged_exposures; }

  unsigned nInt() const { return m_tint.size(); }
  unsigned nEInt() const { return m_eint.size(); }
  // Intervals with: start,stop=MET-t0(), x1=dLivetime, x2=sum(dLivetime)
  const std::vector<TimeInterval>& tInt() const { return m_tint; }
  // Intervals with: start,stop=exposure, x1=aeff, x2=dExposure
  const std::vector<ExposureInterval>& eInt(unsigned ie) const { return m_eint[ie]; }
  unsigned findTInterval(double t) const;
  unsigned findEInterval(double t, unsigned ie) const;
  const std::vector<std::pair<double,double> >& getExpLC() const { return m_lcexp; }

  double t0() const { return m_t0; }

  bool partiallySerializeToBlob(BLOBSerializer& s) const;
  bool partiallyUnserializeFromBlob(BLOBUnserializer& s);

private:
  bool                                          m_t0set;
  double                                        m_t0;
  std::vector<TimeInterval>                     m_tint;
  std::vector<std::vector<ExposureInterval> >   m_eint;
  const std::vector<PLIntegratedEA*>&           m_ea;
  Accumulator                                   m_tacc;
  std::vector<Accumulator>                      m_eacc;
  bool                                          m_merged_exposures;
  double                                        m_lcperiod;
  unsigned                                      m_lcbin;
  Accumulator                                   m_lcacc;
  std::vector<std::pair<double,double> >        m_lcexp;
};

#endif // #ifndef FT2EXP_HPP
