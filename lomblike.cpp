//-*-mode:c++; mode:font-lock;-*-

// lomelike.cpp - Stephen Fegan 
//              - sfegan@llr.in2p3.fr
//              - 2010-07-07
//
// Apply "Lomb"-like Likelihood test for oscillatory model
//
// $Id: lomblike.cpp 3443 2011-12-06 11:56:33Z sfegan $

//#define DUMP_FNC_EVAL

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <irfInterface/IrfsFactory.h>
#include <irfLoader/Loader.h>

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"
#include "VERITAS/VSFileUtility.hpp"
#include "VERITAS/VSDataConverter.hpp"

#include "FT1ROIFilter.hpp"
#include "FT2Exp.hpp"
#include "PLIntegratedEA.hpp"
#include "GTI.hpp"
#include "FITS.hpp"
#include "Util.hpp"
#include "MyMinuit.hpp"
#include "Accumulator.hpp"
#include "MyMinuit.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

// ****************************************************************************
//
// HARMINIC SERIES LIKELIHOOD
//
// ****************************************************************************

class HSLike: public Optimizable
{
public:
  class HarmonicSeries
  {
  public:
    enum ModelComponent { MC_SRC, MC_BKG };

    HarmonicSeries():
      nharmonic(1), freq(0), amp(nharmonic,0), phase(nharmonic,0), 
      amp_free(false), phase_free(false), model_cpt(MC_SRC) { }
    HarmonicSeries(double _amp0, ModelComponent _model_cpt = MC_SRC):
      nharmonic(1), freq(0), amp(nharmonic,_amp0), phase(nharmonic,0), 
      amp_free(true), phase_free(false), model_cpt(_model_cpt) { }
    HarmonicSeries(unsigned _nharmonic, double _freq, 
		   ModelComponent _model_cpt = MC_SRC,
		   double _amp0 = 0, double _phase0 = 0, 
		   bool _amp_free=true, bool _phase_free=true):
      nharmonic(_nharmonic<1?1:_nharmonic), freq(_freq), amp(nharmonic), 
      phase(nharmonic), amp_free(_amp_free), phase_free(_phase_free),
      model_cpt(_model_cpt)
    { 
      amp[0]   = _amp0; phase[0] = _phase0;
      if(freq==0.0) { phase_free=false; nharmonic=1; 
	amp.resize(1); phase.resize(1); phase[0] = 0; }
    }
    
    unsigned               nharmonic;
    double                 freq;
    std::vector<double>    amp;
    std::vector<double>    phase;         // 0 to 1
    bool                   amp_free;
    bool                   phase_free;
    ModelComponent         model_cpt;

    double angfreq(unsigned ih=0) const { return 2.0*M_PI*freq*(ih+1.0); }
    double period(unsigned ih=0) const { return 1.0/freq/(ih+1.0); }
    double angphase(unsigned ih=0) const { return 2.0*M_PI*phase.at(ih); }
    double tau(unsigned ih=0) const { return phase.at(ih)/freq/(ih+1.0); }
  };

  class IntegralFunction
  {
  public:
    IntegralFunction(): 
      clocktime(), livetime(), exposure(), npred(), nmeas(),
      npred_src(), nmeas_prob() {}
    IntegralFunction(double _clocktime, double _livetime, double _exposure,
		     double _npred, double _nmeas, 
		     double _npred_src, double _nmeas_prob): 
      clocktime(_clocktime), livetime(_livetime), exposure(_exposure),
      npred(_npred), nmeas(_nmeas), npred_src(_npred_src), 
      nmeas_prob(_nmeas_prob) {}
    double clocktime;
    double livetime;
    double exposure;
    double npred;
    double nmeas;
    double npred_src;
    double nmeas_prob;
  };

  enum Mode { HMM_SRC_ONLY, HMM_SRCBKG_PROB };

  HSLike(const FT2Exp& ft2exp, const std::vector<FT1AEff>& events,
	 Mode hmm = HMM_SRC_ONLY);
  virtual ~HSLike();

  unsigned addConstant(double amp = 1.0, 
		       HarmonicSeries::ModelComponent mc = 
		       HarmonicSeries::MC_SRC);
  unsigned addHarmonicSeries(unsigned nharmonic, double freq,
			     HarmonicSeries::ModelComponent mc = 
			     HarmonicSeries::MC_SRC,
			     double amp0 = 0, double phi0 = 0);
  unsigned addHarmonicSeries(const HarmonicSeries& hs);
  void replaceHarmonicSeries(unsigned ihs,
			     unsigned nharmonic, double freq,
			     HarmonicSeries::ModelComponent mc = 
			     HarmonicSeries::MC_SRC,
			     double amp0 = 0, double phi0 = 0);
  void replaceHarmonicSeries(unsigned ihs, const HarmonicSeries& hs);
  void deleteHarmonicSeries(unsigned ihs);

  unsigned nHarmonicSeries() const { return m_hsdata.size(); }

  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& params);
  virtual double dfdp(const std::vector<double>& params, unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);

  double srcNorm() const { return m_src_norm; }
  double bkgNorm() const { return m_bkg_norm; }
  double norm(HarmonicSeries::ModelComponent mc) const 
  {
    return (mc == HarmonicSeries::MC_SRC)?srcNorm():bkgNorm();
  }

  bool setAB(unsigned ihs, unsigned ih, double a, double b);
  void enableHS(unsigned ihs);
  void disableHS(unsigned ihs);

  bool setP(const std::vector<double>& params);
  std::vector<double> getP() const;

  void getIntegralFunction(std::vector<IntegralFunction>& fn) const;

private:

  class HSData
  {
  public:
    typedef HarmonicSeries::ModelComponent ModelComponent;

    class HData
    {
    public:
      double              phi0;
      double              cphi0;
      double              sphi0;
      double              a;
      double              b;
      double              int_ecwt;
      double              int_eswt;
      std::vector<double> cwti;
      std::vector<double> swti;
    };
    
    HSData(const HarmonicSeries& hs, unsigned nevent, Mode hmm);
    unsigned nFreeParam() const 
    {
      if(!enabled)return 0;
      unsigned nfree = 0;
      if(amp_free)nfree += nharmonic;
      if(phase_free)nfree += nharmonic;
      return nfree;
    }
    std::string pName(unsigned ihsdata, unsigned iparam);
    double plo(unsigned iparam);
    double phi(unsigned iparam);
    double pval(unsigned iparam);
    bool setPval(unsigned iparam, double x);
    double getPval(unsigned iparam) const;
    void updateData(double int_exposure, const FT2Exp& ft2exp,
		    const std::vector<double>& evtime, bool force_full_update);
    void accumulateMWNPF(unsigned i0, unsigned iN, const FT2Exp& ft2exp,
			 std::vector<Accumulator>& ac, 
			 std::vector<Accumulator>& as);
    
    bool                  enabled;
    ModelComponent        model_cpt;
    unsigned              nharmonic;
    bool                  amp_free;
    bool                  phase_free;
    double                angfreq;
    std::vector<HData>    harmonic;
  };

  inline double dfda(HSData& hsd, unsigned ih) const;
  inline double dfdb(HSData& hsd, unsigned ih) const;
  inline double getEvPr(unsigned ievent, 
			Accumulator& a_src, Accumulator& a_bkg) const;
  inline double getEvPr(unsigned ievent) const
  { Accumulator a_src, a_bkg; return getEvPr(ievent,a_src,a_bkg); }
  inline void nPredFactors(double &src_factor, double& bkg_factor) const;
  inline double nPred() const;
  void recalcEvPr();

  template<typename T> static 
  T* chooseSB(HarmonicSeries::ModelComponent mc, T* s, T* b)
  { return mc==HarmonicSeries::MC_SRC?s:b; }
  template<typename T> static 
  T* chooseSB(HarmonicSeries::ModelComponent mc, T& s, T& b)
  { return chooseSB(mc,&s,&b); }
  template<typename T> static 
  T* chooseSB(const HSData& hsd, T* s, T*b)
  { return chooseSB(hsd.model_cpt,s,b); }
  template<typename T> static 
  T* chooseSB(const HSData& hsd, T& s, T&b)
  { return chooseSB(hsd.model_cpt,s,b); }
  
  Mode                          m_hmm;
  Mode                          m_hmm_req;
  const FT2Exp&                 m_ft2exp;
  const std::vector<FT1AEff>&   m_events;

  double                        m_exposure;
  double                        m_src_norm;
  double                        m_bkg_norm;

  std::vector<HSData*>          m_hsdata;
  std::vector<double>           m_evtime;
  //  std::vector<double>           m_ev_src_weight;
  //  std::vector<double>           m_ev_bkg_weight;
  std::vector<double>           m_ev_spatial_prob0;
  std::vector<double>           m_ev_density;
		      
  double                        m_src_npred_factor;
  double                        m_bkg_npred_factor;
};

HSLike::HSData::
HSData(const HSLike::HarmonicSeries& hs, unsigned nevent,
       HSLike::Mode hmm):
  enabled(true), model_cpt(hs.model_cpt), nharmonic(hs.nharmonic),
  amp_free(hs.amp_free), phase_free(hs.amp_free & hs.phase_free),
  angfreq(hs.angfreq()), harmonic(hs.nharmonic)
{
  if(hmm == HMM_SRC_ONLY)model_cpt = HarmonicSeries::MC_SRC;
  for(unsigned ih=0;ih<hs.nharmonic;ih++)
    {
      HData& h(harmonic[ih]);
      h.phi0        = hs.angphase(ih);
      h.cphi0       = std::cos(h.phi0);
      h.sphi0       = std::sin(h.phi0);
      h.a           = hs.amp[ih];
      h.b           = 0;
      h.int_ecwt    = 1;
      h.int_eswt    = 0;
      h.cwti.resize(nevent,1);
      h.swti.resize(nevent,0);
    }
}

std::string HSLike::HSData::pName(unsigned ihsdata, unsigned iparam)
{
  std::string name;
  if(nharmonic>1)
    {
      unsigned ih = iparam;
      if(phase_free)
	{
	  ih/=2;
	  if(iparam%2 == 0)name = std::string("a");
	  else name = std::string("b");
	}
      else name = std::string("a");
      name += VSDataConverter::toString(ihsdata);
      name += std::string("x");
      name += VSDataConverter::toString(ih+1);
    }
  else
    {
      if(iparam==0)
	name = std::string("a")+VSDataConverter::toString(ihsdata);
      else
	name = std::string("b")+VSDataConverter::toString(ihsdata);
    }
  return name;
}

double HSLike::HSData::plo(unsigned iparam)
{
  if(angfreq==0)return 0;
  return -2.0;
}

double HSLike::HSData::phi(unsigned iparam)
{
  return 2.0;
}

double HSLike::HSData::pval(unsigned iparam)
{
  unsigned ih = iparam;
  if(phase_free)
    {
      ih/=2;
      if(iparam%2 == 0)return harmonic[ih].a;
      else return harmonic[ih].b;
    }
  else
    {
      return harmonic[ih].a;
    }
}

bool HSLike::HSData::setPval(unsigned iparam, double x)
{
  bool changed(true);
  unsigned ih = iparam;
  if(phase_free)
    {
      ih/=2;
      if(iparam%2 == 0){ changed = (harmonic[ih].a != x); harmonic[ih].a = x; }
      else { changed = (harmonic[ih].b != x); harmonic[ih].b = x; }
    }
  else { changed = (harmonic[ih].a != x); harmonic[ih].a = x; }
  return changed;
}

double HSLike::HSData::getPval(unsigned iparam) const
{
  unsigned ih = iparam;
  if(phase_free)
    {
      ih/=2;
      if(iparam%2 == 0)return harmonic[ih].a;
      else return harmonic[ih].b;
    }
  else return harmonic[ih].a;
  assert(0);
}

void HSLike::HSData::
accumulateMWNPF(unsigned i0, unsigned iN, const FT2Exp& ft2exp,
		std::vector<Accumulator>& ac, std::vector<Accumulator>& as)
{
  // Accumulate model weighted NPred factor for this harmonic series
  const unsigned nh = harmonic.size();
  const double wi = angfreq;
  const unsigned ne = ft2exp.nEInt();
  if(wi == 0)
    {
      for(unsigned ih=0;ih<nh;ih++)
	{
	  double ej = 0;
	  for(unsigned ie=0;ie<ne;ie++)
	    if(!ft2exp.eInt(ie).empty())
	      {
		if(iN>0)ej += ft2exp.eInt(ie)[iN-1].exp_int_stop;
		if(i0>0)ej -= ft2exp.eInt(ie)[i0].exp_int_start;
	      }
	  ac[ih].add(ej);
	  //as[if].add(0);
	}
    }
  else
    {
      for(unsigned jint=i0;jint<iN;jint++)
	{
	  const double witj = 
	    wi*0.5*(ft2exp.tInt()[jint].clock_time_start +
		    ft2exp.tInt()[jint].clock_time_stop);
	  double ej = 0;
	  for(unsigned ie=0;ie<ft2exp.nEInt();ie++)
	    if(!ft2exp.eInt(ie).empty())
	      ej += ft2exp.eInt(ie)[jint].dexp;

	  const double cwitj = std::cos(witj);
	  const double switj = std::sin(witj);
	  ac.front().add(ej*cwitj);
	  as.front().add(ej*switj);
	  if(nh>1)
	    {
	      double ckwitj = cwitj;
	      double skwitj = switj;
	      for(unsigned kh=1;kh<nh;kh++)
		{
		  const double _ckwitj = ckwitj*cwitj - skwitj*switj;
		  const double _skwitj = skwitj*cwitj + ckwitj*switj;
		  ckwitj = _ckwitj;
		  skwitj = _skwitj;
		  ac[kh].add(ej*ckwitj);
		  as[kh].add(ej*skwitj);
		}
	    }
	}
    }
}
				   
void HSLike::HSData::
updateData(double int_exposure, const FT2Exp& ft2exp,
	   const std::vector<double>& evtime, bool force_full_update)
{
  const unsigned nh = harmonic.size();
  if(force_full_update)
    {
      const double wi = angfreq;

      std::vector<Accumulator> ac(nh);
      std::vector<Accumulator> as(nh);
      accumulateMWNPF(0,ft2exp.nInt(),ft2exp,ac,as);
      for(unsigned ih=0;ih<nh;ih++)
	{
	  HData& hd(harmonic[ih]);
	  hd.int_ecwt = 
	    (ac[ih].sum()*hd.cphi0 + as[ih].sum()*hd.sphi0)/int_exposure;
	  hd.int_eswt = 
	    (as[ih].sum()*hd.cphi0 - ac[ih].sum()*hd.sphi0)/int_exposure;
	  //	  std::cout << "ECWT: " << wi << ' ' 
	  //		    << hd.int_ecwt << ' ' << hd.int_eswt << '\n';
	}

      if(wi == 0)
	{
	  for(unsigned jevent=0;jevent<evtime.size();jevent++)
	    {
	      for(unsigned ih=0;ih<nh;ih++)
		{
		  HData& hd(harmonic[ih]);
		  hd.cwti[jevent] = 1.0;
		  hd.swti[jevent] = 0.0;
		}
	    }
	}
      else
	{
	  for(unsigned jevent=0;jevent<evtime.size();jevent++)
	    {
	      const double witj = wi*evtime[jevent];
	      const double cwitj = std::cos(witj);
	      const double switj = std::sin(witj);
	      HData& hd0(harmonic.front());
	      hd0.cwti[jevent] = cwitj*hd0.cphi0 + switj*hd0.sphi0;
	      hd0.swti[jevent] = switj*hd0.cphi0 - cwitj*hd0.sphi0;
	      if(nh>1)
		{
		  double ckwitj = cwitj;
		  double skwitj = switj;
		  for(unsigned kh=1;kh<nh;kh++)
		    {
		      HData& hd(harmonic[kh]);
		      const double _ckwitj = ckwitj*cwitj - skwitj*switj;
		      const double _skwitj = skwitj*cwitj + ckwitj*switj;
		      ckwitj = _ckwitj;
		      skwitj = _skwitj;
		      hd.cwti[jevent] = ckwitj*hd.cphi0 + skwitj*hd.sphi0;
		      hd.swti[jevent] = skwitj*hd.cphi0 - ckwitj*hd.sphi0;
		    }
		}
	    }
	}
    }
}

HSLike::
HSLike(const FT2Exp& ft2exp, const std::vector<FT1AEff>& events,
       Mode hmm):
  Optimizable(), 
  m_hmm(hmm), m_ft2exp(ft2exp), m_events(events), 
  m_exposure(), m_src_norm(), m_bkg_norm(), 
  m_hsdata(), m_evtime(events.size()),
  m_ev_spatial_prob0(), 
  m_ev_density(events.size()), m_src_npred_factor(), m_bkg_npred_factor()
{  
  for(unsigned ie=0;ie<m_ft2exp.nEInt();ie++)
    m_exposure += m_ft2exp.eInt(ie).back().stop; 
  m_src_norm = double(m_events.size())/m_exposure;
  m_bkg_norm = 0;
  for(unsigned ievent=0;ievent<events.size();ievent++)
    m_evtime[ievent] = events[ievent].time-ft2exp.t0();
  if(hmm != HMM_SRC_ONLY)
    {
      Accumulator a;
      m_ev_spatial_prob0.resize(events.size());
      for(unsigned ievent=0;ievent<events.size();ievent++)
	{
	  const double p0 =  events[ievent].x_dbl.at(0);
	  m_ev_spatial_prob0[ievent] = p0;
	  a.add(p0);
	}
      m_src_norm = a.sum()/m_exposure;
      m_bkg_norm = (double(m_events.size())-a.sum())/m_exposure;
    }
}

HSLike::~HSLike()
{
  for(std::vector<HSData*>::iterator ihsd = m_hsdata.begin();
      ihsd != m_hsdata.end(); ihsd++)delete *ihsd;
}

unsigned HSLike::addConstant(double amp,
			     HarmonicSeries::ModelComponent mc)
{
  return addHarmonicSeries(HarmonicSeries(amp,mc));
}

unsigned HSLike::addHarmonicSeries(unsigned nharmonic, double freq,
				   HarmonicSeries::ModelComponent mc,
				   double amp0, double phi0)
{
  return addHarmonicSeries(HarmonicSeries(nharmonic,freq,mc,amp0,phi0));
}

unsigned HSLike::addHarmonicSeries(const HarmonicSeries& hs)
{
  unsigned ihs = m_hsdata.size();
  m_hsdata.push_back(new HSData(hs, m_events.size(), m_hmm));  
  m_hsdata[ihs]->updateData(m_exposure, m_ft2exp, m_evtime, true);
  recalcEvPr();						
  return ihs;
}

void HSLike::replaceHarmonicSeries(unsigned ihs,
				   unsigned nharmonic, double freq,
				   HarmonicSeries::ModelComponent mc,
				   double amp0, double phi0)
{
  replaceHarmonicSeries(ihs,HarmonicSeries(nharmonic,freq,mc,amp0,phi0));
}

void HSLike::replaceHarmonicSeries(unsigned ihs, const HarmonicSeries& hs)
{
  assert(ihs<m_hsdata.size());
  delete m_hsdata[ihs];
  m_hsdata[ihs] = new HSData(hs, m_events.size(), m_hmm);
  m_hsdata[ihs]->updateData(m_exposure, m_ft2exp, m_evtime, true);
  recalcEvPr();
}

void HSLike::deleteHarmonicSeries(unsigned ihs)
{
  assert(ihs<m_hsdata.size());
  delete m_hsdata[ihs];
  m_hsdata.erase(m_hsdata.begin()+ihs);
  recalcEvPr();
}

bool HSLike::setAB(unsigned ihs, unsigned ih, double a, double b)
{
  assert(ihs<m_hsdata.size());
  HSData& hsd(*m_hsdata[ihs]);
  assert(ih<hsd.harmonic.size());
  assert((hsd.angfreq!=0.0&&!hsd.phase_free)||b==0.0); // WTF?
  HSData::HData& hd(hsd.harmonic[ih]);
  bool changed = false;
  if(hd.a != a){ hd.a = a; changed |= true; }
  if(hd.b != b){ hd.b = b; changed |= true; }
  if(changed)
    {
      hsd.updateData(m_exposure,m_ft2exp,m_evtime,false);
      recalcEvPr();
    }
  return changed;
}

void HSLike::enableHS(unsigned ihs)
{
  assert(ihs<m_hsdata.size());
  if(!m_hsdata[ihs]->enabled)
    {
      m_hsdata[ihs]->enabled=true;
      recalcEvPr();
    }
}

void HSLike::disableHS(unsigned ihs)
{
  assert(ihs<m_hsdata.size());
  if(m_hsdata[ihs]->enabled)
    {
      m_hsdata[ihs]->enabled=false;
      recalcEvPr();
    }
}

void writeP(const std::vector<double>& params)
{
  for(unsigned ip=0;ip<params.size();ip++)
    {
      if(ip!=0)std::cout << ' ';
      std::cout << params[ip];
    }
}

bool HSLike::setP(const std::vector<double>& params)
{
  bool changed = false;
  const unsigned nhs = m_hsdata.size();
  unsigned ifree = 0;
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      bool hs_changed = false;
      HSData& hsd(*m_hsdata[ihs]);
      unsigned nfp = hsd.nFreeParam();
      for(unsigned ifp=0;ifp<nfp;ifp++,ifree++)
	hs_changed |= hsd.setPval(ifp,params[ifree]);
      if(hs_changed)hsd.updateData(m_exposure,m_ft2exp,m_evtime,false);
      changed |= hs_changed;
    }
  if(changed)recalcEvPr();
  return changed;
}

std::vector<double> HSLike::getP() const
{
  std::vector<double> p;
  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      unsigned nfp = hsd.nFreeParam();
      for(unsigned ifp=0;ifp<nfp;ifp++)p.push_back(hsd.getPval(ifp));
    }
  return p;
}

void HSLike::getIntegralFunction(std::vector<IntegralFunction>& fn) const
{
  fn.clear();
  fn.push_back(IntegralFunction());
  const unsigned nevent = m_events.size();
  const unsigned ne = m_ft2exp.nEInt();
  const unsigned nhs = m_hsdata.size();
  unsigned iint=0;
  Accumulator a_npred_src;
  Accumulator a_npred_bkg;
  Accumulator a_nmeas_prob;
  for(unsigned ievent=0;ievent<nevent;ievent++)
    {
      if(m_hmm == HMM_SRC_ONLY)a_nmeas_prob.add(1.0);      
      else 
	{
	  Accumulator a_src, a_bkg;
	  const double rho = getEvPr(ievent, a_src, a_bkg);
	  a_nmeas_prob.add(m_ev_spatial_prob0[ievent]*a_src.sum()/rho);
	}

      const double ev_time = m_events[ievent].time-m_ft2exp.t0();
      const unsigned jint = m_ft2exp.findTInterval(ev_time);

      double exposure = 0;
      for(unsigned ie=0;ie<ne;ie++)
	exposure += m_ft2exp.eInt(ie)[jint].exp_int_stop;

      for(unsigned ihs=0;ihs<nhs;ihs++)
	{
	  HSData& hsd(*m_hsdata[ihs]);
	  if(!hsd.enabled)continue;
	  Accumulator& a(*chooseSB(hsd,a_npred_src,a_npred_bkg));
	  const double hsdnorm = norm(hsd.model_cpt);

	  const unsigned nh(hsd.harmonic.size());
	  std::vector<Accumulator> ac(nh);
	  std::vector<Accumulator> as(nh);
	  hsd.accumulateMWNPF(iint,jint+1,m_ft2exp,ac,as);
	  for(unsigned ih=0;ih<nh;ih++)
	    {
	      HSData::HData& hd(hsd.harmonic[ih]);
	      a.add(hsdnorm*(hd.a*ac[ih].sum() + hd.b*as[ih].sum()));
	    }
	}
      IntegralFunction ifn;
      ifn.clocktime     = m_ft2exp.tInt()[jint].clock_time_stop;
      ifn.livetime      = m_ft2exp.tInt()[jint].livetime_int_stop;
      ifn.exposure      = exposure;
      ifn.npred         = a_npred_src.sum()+a_npred_bkg.sum();
      ifn.nmeas         = ievent+1;
      ifn.npred_src     = a_npred_src.sum();
      ifn.nmeas_prob    = a_nmeas_prob.sum();
      fn.push_back(ifn);      
      iint = jint+1;
    }

  fn.push_back(IntegralFunction(m_ft2exp.tInt().back().clock_time_stop,
				m_ft2exp.tInt().back().livetime_int_stop,
				m_exposure, nPred(), nevent,
				m_exposure*m_src_npred_factor*m_src_norm,
				a_nmeas_prob.sum()));
}

void HSLike::recalcEvPr()
{
  nPredFactors(m_src_npred_factor,m_bkg_npred_factor);

  unsigned nevent = m_ev_density.size();
  bool warned = false;
  for(unsigned ievent=0;ievent<nevent;ievent++)
    {
      const double rho = getEvPr(ievent);
      m_ev_density[ievent] = rho<0 ? 0 : rho;
      if(rho<0 && warned==false)
	{
	  std::cout << "Warning: rho=" << rho << '\n';
	  warned = true;
	}
    }
}

unsigned HSLike::numParam()
{
  unsigned nfree = 0;
  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      nfree += hsd.nFreeParam();
    }
  return nfree;
}

std::string HSLike::pName(unsigned iparam)
{
  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      unsigned nfree = hsd.nFreeParam();
      if(iparam>=nfree)iparam-=nfree;
      else return hsd.pName(ihs,iparam);
    }
  assert(0);
}
 
bool HSLike::canCalcDFDP()
{
#if 1
  return true;
#else
  return false;
#endif
}
 
inline double HSLike::getEvPr(unsigned ievent,
			      Accumulator& a_src, Accumulator& a_bkg) const
{
  //Accumulator a_src;
  //Accumulator a_bkg;
  const unsigned nhs(m_hsdata.size());
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      if(!hsd.enabled)continue;
      Accumulator* a(chooseSB(hsd,a_src,a_bkg));

      // If statement not strictly necessary - 2nd formula would do for both
      if(hsd.angfreq == 0.0)
	{
	  a->add(hsd.harmonic.front().a);
	}
      else
	{
	  const unsigned nh(hsd.harmonic.size());
	  for(unsigned ih=0;ih<nh;ih++)
	    {
	      HSData::HData& hd(hsd.harmonic[ih]);
	      a->add(hd.a*hd.cwti[ievent] + hd.b*hd.swti[ievent]);
	    }
	}
    }

  double rho_event = 0;
  switch(m_hmm)
    {
    case HMM_SRC_ONLY:
      rho_event = a_src.sum();
      break;
    case HMM_SRCBKG_PROB:
      {
	const double p0 = m_ev_spatial_prob0[ievent];
	rho_event = (1.0-p0)*a_bkg.sum() + p0*a_src.sum();
      }
      break;
    }
  return rho_event;
}

inline void HSLike::nPredFactors(double& src_factor, double& bkg_factor) const
{
  Accumulator a_src;
  Accumulator a_bkg;
  const unsigned nhs(m_hsdata.size());
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      if(!hsd.enabled)continue;

      Accumulator* a(chooseSB(hsd,a_src,a_bkg));

      // If statement not strictly necessary - 2nd formula would do for both
      if(hsd.angfreq == 0.0)
	{
	  a->add(hsd.harmonic.front().a);
	}
      else
	{
	  const unsigned nh(hsd.harmonic.size());
	  for(unsigned ih=0;ih<nh;ih++)
	    {
	      HSData::HData& hd(hsd.harmonic[ih]);
	      a->add(hd.a*hd.int_ecwt + hd.b*hd.int_eswt);
	    }
	}
    }

  src_factor = a_src.sum();
  bkg_factor = a_bkg.sum();
}

inline double HSLike::nPred() const
{
  return m_exposure*(m_src_npred_factor*m_src_norm + 
		     m_bkg_npred_factor*m_bkg_norm);
}

double HSLike::f(const std::vector<double>& params)
{
  setP(params);
  Accumulator a;

  // Data sum
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      const double rho = m_ev_density[ievent];
      a.add(std::log(rho));
    }
  
  // Model sum
  const double npred = nPred();

#ifdef DUMP_FNC_EVAL
  Accumulator px_src;
  Accumulator px_bkg;
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      const double p = m_ev_spatial_prob0[ievent];
      px_src.add(p);
      px_bkg.add(1.0-p);
    }
  writeP(params); 
  std::cout << " F: " << -(a.sum() - npred) << ' ' 
	    << px_src.sum() << ' ' 
	    << m_src_npred_factor*m_src_norm*m_exposure << ' '
	    << px_bkg.sum() << ' ' 
	    << m_bkg_npred_factor*m_bkg_norm*m_exposure << '\n';
#endif
  return -(a.sum() - npred);
}

double HSLike::dfda(HSData& hsd, unsigned ih) const
{
  Accumulator a;
 
  // Data sum
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      const double rho_event = m_ev_density[ievent];
      double devsum_da = 1.0/rho_event;
      if(hsd.angfreq != 0.0)
	devsum_da *= hsd.harmonic[ih].cwti[ievent];

      if(hsd.model_cpt == HarmonicSeries::MC_BKG)
	devsum_da *= 1.0 - m_ev_spatial_prob0[ievent];
      else if(m_hmm == HMM_SRCBKG_PROB)
	devsum_da *= m_ev_spatial_prob0[ievent];
      a.add(devsum_da);
    }

  // Model sum
  const double hsdnorm = norm(hsd.model_cpt);
  double dnpred_da(0);
  if(hsd.angfreq == 0.0)dnpred_da = hsdnorm*m_exposure;
  else dnpred_da = hsdnorm*hsd.harmonic[ih].int_ecwt*m_exposure;

#ifdef DUMP_FNC_EVAL
  std::cout << " dFdA: " << -(a.sum() - dnpred_da) << '\n';
#endif
  return -(a.sum() - dnpred_da);
}

double HSLike::dfdb(HSData& hsd, unsigned ih) const
{
  Accumulator a;
 
  // Temporary
  assert(hsd.angfreq != 0.0 && hsd.phase_free);

  // Data sum
  for(unsigned ievent=0;ievent<m_events.size();ievent++)
    {
      const double rho_event = m_ev_density[ievent];
      double devsum_db = hsd.harmonic[ih].swti[ievent]/rho_event;

      if(hsd.model_cpt == HarmonicSeries::MC_BKG)
	devsum_db *= 1.0 - m_ev_spatial_prob0[ievent];
      else if(m_hmm == HMM_SRCBKG_PROB)
	devsum_db *= m_ev_spatial_prob0[ievent];
      a.add(devsum_db);
    }

  // Model sum
  const double hsdnorm = norm(hsd.model_cpt);
  double dnpred_db = hsdnorm*hsd.harmonic[ih].int_eswt*m_exposure;

#ifdef DUMP_FNC_EVAL
  std::cout << " dFdB: " << -(a.sum() - dnpred_db) << '\n';
#endif
  return -(a.sum() - dnpred_db);
}

double HSLike::dfdp(const std::vector<double>& params, unsigned iparam)
{
  setP(params);
#ifdef DUMP_FNC_EVAL
  writeP(params); 
#endif

  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      if(!hsd.enabled)continue;
      unsigned nfree = hsd.nFreeParam();
      if(iparam>=nfree)iparam -= nfree;
      else 
	{
	  if(hsd.phase_free)
	    {
	      unsigned ih = iparam/2;
	      if(iparam%2 == 0)return dfda(hsd,ih);
	      else return dfdb(hsd,ih);
	    }
	  else return dfda(hsd,iparam);
	}
    }
  assert(0);
}

double HSLike::p0(unsigned iparam)
{
  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      if(!hsd.enabled)continue;
      unsigned nfree = hsd.nFreeParam();
      if(iparam>=nfree)iparam-=nfree;
      else return hsd.pval(iparam);
    }
  assert(0);
}

double HSLike::plo(unsigned iparam)
{
  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      if(!hsd.enabled)continue;
      unsigned nfree = hsd.nFreeParam();
      if(iparam>=nfree)iparam-=nfree;
      else return hsd.plo(iparam);
    }
  assert(0);
}

double HSLike::phi(unsigned iparam)
{
  const unsigned nhs = m_hsdata.size();
  for(unsigned ihs=0;ihs<nhs;ihs++)
    {
      HSData& hsd(*m_hsdata[ihs]);
      if(!hsd.enabled)continue;
      unsigned nfree = hsd.nFreeParam();
      if(iparam>=nfree)iparam-=nfree;
      else return hsd.phi(iparam);
    }
  assert(0);
}

class Scanner
{
public:
  class Scan
  {

  };

  class KnownFrequency
  {
  public:
    KnownFrequency(): 
      freq(0), nharmonic(1), cpt(HSLike::HarmonicSeries::MC_SRC) { }
    KnownFrequency(double _freq, unsigned _nharmonic = 1,
		   HSLike::HarmonicSeries::ModelComponent _cpt = 
		   HSLike::HarmonicSeries::MC_SRC):
      freq(_freq), nharmonic(_nharmonic), cpt(_cpt) { }
    double freq;
    unsigned nharmonic;
    HSLike::HarmonicSeries::ModelComponent cpt;
  };

  Scanner(double scan0, double scanN, unsigned nscan, unsigned nharmonic,
	  const FT2Exp& ft2exp, unsigned verbose, bool log_scan = false,
	  HSLike::Mode mode = HSLike::HMM_SRC_ONLY):
    m_scan0(scan0), m_scanN(scanN), m_nscan(nscan), m_nharmonic(nharmonic),
    m_ft2exp(ft2exp), m_verbose(verbose), m_log_scan(log_scan), m_mode(mode),
    m_null_logl(), m_null_s(), m_null_ds(), m_null_b(), m_null_db() { }
    
  void scan(const std::vector<FT1AEff>& events,
	    const std::vector<KnownFrequency>& known_freq
	    = std::vector<KnownFrequency>());
  
private:
  double         m_scan0;
  double         m_scanN;
  unsigned       m_nscan;
  unsigned       m_nharmonic;
  const FT2Exp&  m_ft2exp;
  unsigned       m_verbose;
  bool           m_log_scan;
  HSLike::Mode   m_mode;

  double         m_null_logl;
  double         m_null_s;
  double         m_null_ds;
  double         m_null_b;
  double         m_null_db;
};

void Scanner::scan(const std::vector<FT1AEff>& events,
		   const std::vector<KnownFrequency>& known_freq)
{
  double exposure(0);
  for(unsigned ie=0;ie<m_ft2exp.nEInt();ie++)
    exposure += m_ft2exp.eInt(ie).back().stop;

  HSLike m(m_ft2exp, events, m_mode);
  if(m_mode != HSLike::HMM_SRC_ONLY)
    m.addConstant(1.0,HSLike::HarmonicSeries::MC_BKG);
  m.addConstant(1.0,HSLike::HarmonicSeries::MC_SRC);

  for(unsigned ikf = 0;ikf<known_freq.size();ikf++)
    m.addHarmonicSeries(known_freq[ikf].nharmonic, known_freq[ikf].freq, 
			known_freq[ikf].cpt);

  std::vector<double> default_param_values;

  if(1)
    {
      MyMinuit minuit(&m,std::max(int(m_verbose)-4,-1));
      minuit.minimize(false,1e-3,0.5,"",true,2);
      //minuit.minimize(true,1e-3,0.5,"",true,2);
      m_null_logl = minuit.fVal();
      if(m_mode != HSLike::HMM_SRC_ONLY)
	{
	  m_null_b  = minuit.pVal(0) * m.bkgNorm();
	  m_null_db = minuit.pErr(0) * m.bkgNorm();
	  m_null_s  = minuit.pVal(1) * m.srcNorm();
	  m_null_ds = minuit.pErr(1) * m.srcNorm();
	}
      else
	{
	  m_null_s  = minuit.pVal(0) * m.srcNorm();
	  m_null_ds = minuit.pErr(0) * m.srcNorm();
	}
      for(unsigned ip=0;ip<m.numParam();ip++)
	default_param_values.push_back(minuit.pVal(ip));
    }

  if(1)
    {
      m.setP(default_param_values);
      std::vector<HSLike::IntegralFunction> ifn;
      m.getIntegralFunction(ifn);
      std::ofstream ifnstream("ifn.dat");
      for(unsigned iifn=0;iifn<ifn.size();iifn++)
	ifnstream 
	  << ifn[iifn].clocktime << ' '
	  << ifn[iifn].livetime << ' '
	  << ifn[iifn].exposure << ' '
	  << ifn[iifn].npred << ' '
	  << ifn[iifn].nmeas << ' '
	  << ifn[iifn].npred_src << ' '
	  << ifn[iifn].nmeas_prob << '\n';
    }

  if(m_verbose>=2)
    {
      if(m_mode != HSLike::HMM_SRC_ONLY)
	std::cout 
	  << "FIT: Const model B=" 
	  << m_null_b << "+/-" << m_null_db << "  S="
	  << m_null_s << "+/-" << m_null_ds << "  LogL=" 
	  << m_null_logl << '\n';
      else
	std::cout 
	  << "FIT: Const model S=" << m_null_s << "+/-" << m_null_ds << ' '
	  << "LogL=" << m_null_logl << '\n';
    }

  const unsigned iscancpt = m.nHarmonicSeries();
  const unsigned iscanpar = default_param_values.size();
  for(unsigned iscan=0;iscan<m_nscan;iscan++)
    {
      double xscan = m_scan0;
      if(m_nscan>1)
	{
	  if(m_log_scan)
	    {
	      const double x0 = -std::log(m_scanN);
	      const double xN = -std::log(m_scan0);
	      xscan = x0+(xN-x0)/double(m_nscan-1)*double(m_nscan-1-iscan);
	      xscan = exp(-xscan);
	    }
	  else
	    {
	      const double x0 = m_scan0;
	      const double xN = m_scanN;
	      xscan = x0+(xN-x0)/double(m_nscan-1)*double(iscan);
	    }
	}

      if(iscan==0)
	{
	  m.addHarmonicSeries(m_nharmonic,xscan,
			      HSLike::HarmonicSeries::MC_SRC,0,0);
	  default_param_values.resize(m.numParam());
	}
      else 
	{
	  m.replaceHarmonicSeries(iscancpt,m_nharmonic,xscan,
				  HSLike::HarmonicSeries::MC_SRC,0,0);
	}
      m.setP(default_param_values);

      //if(m_phase_scan)harmonic[1].phase = xscan;
      //else harmonic[1].freq = xscan;

      MyMinuit minuit(&m,std::max(int(m_verbose)-4,-1));

      try
	{
	  minuit.minimize(false,1e-3,0.5,"",true,2);
	}
      catch(const std::string& x)
	{
	  std::cout << x;
	}
      
      if(m_verbose>=1)
	{
	  const double a  = minuit.pVal(iscanpar);
	  const double da = minuit.pErr(iscanpar);
	  const double b  = minuit.pVal(iscanpar+1);
	  const double db = minuit.pErr(iscanpar+1);
	  
	  std::cout << xscan*86400.0 << ' ' 
		    << 1.0/(86400.0*xscan) << ' ';
	  std::cout 
	    << m_null_logl-minuit.fVal() << ' '
	    << a*m.srcNorm() << ' ' 
	    << da*m.srcNorm() << ' '
	    << b*m.srcNorm() << ' ' 
	    << db*m.srcNorm() << ' '
	    << std::sqrt(a*a + b*b)*m.srcNorm() << ' '
	    << std::atan2(b,a)/2.0/M_PI << '\n';
	}
    }
}

// ****************************************************************************
//
// MAIN
//
// ****************************************************************************

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] ft1_file ft2_file irfs"
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);

  std::string command_line(*argv);
  for(int iarg=1;iarg<argc;iarg++)
    {
      command_line+=" ";
      command_line+=argv[iarg];
    }

  // --------------------------------------------------------------------------
  // PROCESS OPTIONS
  // --------------------------------------------------------------------------

  VSOptions options(argc, argv, true);
  
  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  int verbose = 1;
  if(options.find("q","Print no messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 0;
  if(options.find("v","Print really verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 2;
  if(options.find("vv","Print extremely verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 3;
  if(options.find("vvv","Print extremely verbose messages during analysis "
		  "and enable verbose output from the minimizer.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 4;
  if(options.find("vvvv","Print extremely verbose messages during analysis "
		  "and enable verbose output from the minimizer.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 5;

  // ---------------------------- SCANNING OPTIONS ----------------------------

  double scan0 = 1.0/100.0;
  double scanN = 1.0/1.0;
  unsigned nscan = 100;;
  unsigned nharmonic = 1;

  bool scan_period = false;
  if(options.find("period","Specify scan bounds in period rather than "
		  "frequency.","scan")
     !=VSOptions::FS_NOT_FOUND)
    scan_period = true;

  bool scan_hours = false;
  if(options.find("hours","Specify scan frequency in 1/hours or period in "
		  "hours (rather than 1/day and day respectively).","scan")
     !=VSOptions::FS_NOT_FOUND)
    scan_hours = true;

  options.findWithValue("scan0", scan0,
			"Lowest frequency [1/day or 1/hour], or period "
			"[day or hour] in scan.","scan");
  options.findWithValue("scanN", scanN,
			"Highest frequency [1/day or 1/hour], or period "
			"[day or hour] in scan.","scan");
  options.findWithValue("nscan", nscan,
			"Number of frequencies (or periods) in scan.","scan");
  options.findWithValue("nharmonic", nharmonic,
			"Number of harmonics of the fundamental frequency "
			"to include in the model when scanning.","scan");

  bool log_scan = false;
  if(options.find("scan_log","Scan frequencies (or phases) logarithmically.",
		  "scan")
     !=VSOptions::FS_NOT_FOUND)
    log_scan = true;
  
  double scan_unit_to_seconds = (scan_hours?3600.0:86400.0);
  if(scan_period)scan0 = 1.0 / scan0, scanN = 1.0 / scanN;
  scan0 /= scan_unit_to_seconds;
  scanN /= scan_unit_to_seconds;
  if(scan0 > scanN) { double _scan0 = scan0; scan0=scanN; scanN=_scan0; }

  std::vector<std::pair<double,unsigned> > skfreq;
  options.findWithValue("sfreq", skfreq,
			"List of known harmonics to include in the "
			"source model. Specify as a comma separated "
			"list of harmonic components: freq/nharmonics. "
			"The frequency should be given in units of [1/day] "
			"(i.e. the \"period\" and \"hours\" options do not "
			"apply to this option).", "scan");  
  std::vector<std::pair<double,unsigned> > bkfreq;
  options.findWithValue("bfreq", bkfreq,
			"List of known harmonics to include in the "
			"source model. Specify as a comma separated "
			"list of harmonic components: freq/nharmonics. "
			"The frequency should be given in units of [1/day] "
			"(i.e. the \"period\" and \"hours\" options do not "
			"apply to this option).", "scan");  

  std::vector<Scanner::KnownFrequency> known_freq;
  for(unsigned ikf=0;ikf<skfreq.size();ikf++)
    known_freq.push_back(Scanner::KnownFrequency(skfreq[ikf].first / 86400.0,
						 skfreq[ikf].second));
  for(unsigned ikf=0;ikf<bkfreq.size();ikf++)
    known_freq.push_back(Scanner::KnownFrequency(bkfreq[ikf].first / 86400.0,
						 bkfreq[ikf].second,
					      HSLike::HarmonicSeries::MC_BKG));

  // ---------------------------- EVENT CUT OPTIONS ---------------------------

  bool has_energy_cuts = false;
  std::pair<double,double> energy_cuts;
  bool has_roi_cuts = false;
  triple<double,double,double> roi_cuts;
  if(options.findWithValue("erange", energy_cuts,
			   "Energy range to consider in analysis (specifies "
			   "cuts on energy in FT1 file). If not provided, "
			   "the full range of energies in the FT1 file is "
			   "used. Specify as emin/emax in MeV.",
			   "cuts") != VSOptions::FS_NOT_FOUND)
    has_energy_cuts = true;
  
  if(options.findWithValue("roi", roi_cuts,
			   "Circular region to consider in analysis "
			   "(specifies cuts on event position in FT1 file). "
			   "If not provided, the full ROI is used. Specify "
			   "as ra/dec/radus in degrees.",
			   "cuts") != VSOptions::FS_NOT_FOUND)
    {
      has_roi_cuts = true;
      roi_cuts.first  = d2r(roi_cuts.first);
      roi_cuts.second = d2r(roi_cuts.second);
      roi_cuts.third  = d2r(roi_cuts.third);
    }

  bool has_psf_roi_cuts = false;
  double psf_roi_cut_frac = 1.0;
  if(options.findWithValue("roi_psf_frac", psf_roi_cut_frac,
			   "PSF-based circular region to consider in analysis "
			   "(specifies cuts on event position in FT1 file). "
			   "If not provided, the full ROI is used. Specify "
			   "the fraction 0<F<1.",
			   "cuts") != VSOptions::FS_NOT_FOUND)
    {
      has_psf_roi_cuts = true;
      if(psf_roi_cut_frac>1)psf_roi_cut_frac=1;
      if(psf_roi_cut_frac<0)psf_roi_cut_frac=0;
    }

  std::string src_prob_name;
  options.findWithValue("src_event_prob", src_prob_name,
			"Use event probabilities for given source. If "
			"selected the src-bkg model is used. Specify the ",
			"name of the source. The source probabilities "
			"for this source, from gtsrcprob, must be in the "
			"FT1 file.", "cuts");

  // -------------------------------- EXPOSURE --------------------------------

  double gamma = 2.0;
  options.findWithValue("gamma", gamma,
			"Spectral index to use in weighting of effective "
			"area.", "exposure");

  bool no_livetime_correction = false;
  if(options.find("no_lt_correction","Do not apply livetime correction (if it "
		  "available in the IRFs)",
		  "exposure")
     !=VSOptions::FS_NOT_FOUND)
    no_livetime_correction = true;

  std::string exp_blob_fn = std::string();
  options.findWithValue("exp_file", exp_blob_fn,
			"Load exposure calculation to a file or load it "
			"back from the file if it exists.","exposure");  

  std::string day_exp_lc_fn = std::string();
  options.findWithValue("day_exp_lc", day_exp_lc_fn,
			"Specify the name of a file in which to write the "
			"lightcurve of daily exposure values. This is mostly "
			"meant as a diagnostics that can be compared with a "
			"curve produced by gtexposure.","exposure");

  // --------------------------------------------------------------------------
  // FINISH OPTIONS PROCESSING
  // --------------------------------------------------------------------------

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
	if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  int arg_req = 3;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string ft1_fn(*argv);
  argc--, argv++;

  std::string ft2_fn(*argv);
  argc--, argv++;

  std::string irf(*argv);
  argc--, argv++;  

  // --------------------------------------------------------------------------
  // GET ROI CUTS
  // --------------------------------------------------------------------------

  FITSHeader ft1head;
  if(verbose>=3)
    std::cout << "Load: FT1 header - " << ft1_fn << '\n';
  ft1head.loadFromFITS(ft1_fn,FT1::tableName());
  
  if(!ft1head.has("NDSKEYS"))
    {
      std::cerr << "Fatal: " << ft1_fn << " does not have ROI cuts\n";
      std::exit(EXIT_FAILURE);
    }

  double roi_emin(0);
  double roi_emax(0);
  double roi_ra(0);
  double roi_dec(0);
  double roi_roi(0);

  unsigned nds(0);
  ft1head.getAs(nds,"NDSKEYS");
  for(unsigned ids=0;ids<nds;ids++)
    {
      std::string type = ft1head.get(ft1head.nk("DSTYP",ids+1));
      if(type == "ENERGY")
	{
	  std::string val = ft1head.get(ft1head.nk("DSVAL",ids+1));
	  std::string::size_type ic = val.find(':');
	  if(ic != val.npos)
	    {
	      VSDataConverter::fromString(roi_emin,val.substr(0,ic));
	      VSDataConverter::fromString(roi_emax,val.substr(ic+1));
	    }
	}
      else if(type == "POS(RA,DEC)")
	{
	  std::string val = ft1head.get(ft1head.nk("DSVAL",ids+1));
	  if(val.substr(0,7)=="CIRCLE(")
	    {
	      std::vector<double> dv;
	      VSDataConverter::fromString(dv,val.substr(7,val.length()-8));
	      roi_ra  = d2r(dv.at(0));
	      roi_dec = d2r(dv.at(1));
	      roi_roi = d2r(dv.at(2));
	    }
	}
    }

  if(verbose>=3)
    {
      std::cout << "FT1: " << roi_emin << " <= E <= " << roi_emax << '\n';
      std::cout << "FT1: |(RA,Dec)-(" 
		<< r2d(roi_ra) << ',' << r2d(roi_dec)
		<< ")| <= " << r2d(roi_roi) << '\n';
    }

  double emin(roi_emin);
  double emax(roi_emax);
  if(has_energy_cuts)
    {
      if(energy_cuts.first<emin)
	{
	  std::cerr << "Minimum energy cut less than FT1 value: " 
		    << energy_cuts.first << " < " << roi_emin << '\n';
	  std::exit(EXIT_FAILURE);
	}

      if(energy_cuts.second>emax)
	{
	  std::cerr << "Maximum energy cut greater than FT1 value: " 
		    << energy_cuts.second << " > " << roi_emax << '\n';
	  std::exit(EXIT_FAILURE);
	}

      emin = energy_cuts.first;
      emax = energy_cuts.second;
    }

  double ra(roi_ra);
  double dec(roi_dec);
  double roi(roi_roi);
  if(has_roi_cuts)
    {
      double d = roi_cuts.third +
	sphere_dist(ra, dec, roi_cuts.first, roi_cuts.second);

      if(d > roi)
	{
	  std::cerr << "Desired ROI is outside FT1 file region: "
		    << r2d(d) << " > " << r2d(roi) << '\n';
	  std::exit(EXIT_FAILURE);
	}
      
      ra  = roi_cuts.first;
      dec = roi_cuts.second;
      roi = roi_cuts.third;
    }

  // --------------------------------------------------------------------------
  // READ IN FT1 GTI
  // --------------------------------------------------------------------------

  GTIRange gti;
  if(verbose>=3)
    std::cout << "Load: GTI - " << ft1_fn << '\n';
  gti.loadGTIsFromFITS(ft1_fn);
  
  if(verbose>=2)
    {
      std::cout << "GTI: Total time " << lrint(gti.totalTime()) 
		<< " seconds in " << gti.nGTIs() << " periods\n";
    }

#if 0
  for(GTIRange::const_iterator igti=gti.begin();igti!=gti.end();igti++)
    std::cout << lrint(igti->t_start) << ' ' << lrint(igti->t_stop) << '\n';
#endif

  // --------------------------------------------------------------------------
  // READ IN IRFS
  // --------------------------------------------------------------------------

  irfLoader::Loader::go();
  const std::vector<std::string> & 
    irf_names(irfLoader::Loader::respIds().find(irf)->second);
  irfInterface::IrfsFactory & factory(*irfInterface::IrfsFactory::instance());
  std::vector<irfInterface::Irfs *> irfs(irf_names.size());
  for (size_t iirf = 0; iirf < irf_names.size(); iirf++)
    irfs[iirf] = factory.create(irf_names[iirf]);

  std::vector<PLIntegratedEA*> eaint(irfs.size());
  for (size_t iirf = 0; iirf < irfs.size(); iirf++)
    eaint[iirf] = new PLIntegratedEA(irfs[iirf], emin, emax, gamma);
  
  // --------------------------------------------------------------------------
  // READ IN FT2
  // --------------------------------------------------------------------------
  
  FT2Exp ft2exp(eaint, gti, ra, dec);
  // ROI CUTS

  bool loaded_from_blob = false;
  if(!exp_blob_fn.empty())
    {
      BLOBUnserializer blob(exp_blob_fn);
      if(verbose>=2)
	std::cout << "Loading exposure BLOB: " << exp_blob_fn;
      if(blob.good())
	loaded_from_blob = ft2exp.partiallyUnserializeFromBlob(blob);
      if(loaded_from_blob && verbose>=2)std::cout << " .. succeeded\n";
      else if(verbose>=2)std::cout << " .. failed\n";
    }

  if(!loaded_from_blob)
    {
      FITSVectorDispatcher<FT2> ft2_dispatcher(&ft2exp);
      if(verbose>=2)
	std::cout << "Load: FT2 - " << ft2_fn << '\n';
      unsigned ft2_ntotal = 
	ft2_dispatcher.dispatchVector(ft2_fn, FT2::tableName());
      if(verbose>=2)
	std::cout
	  << "FT2: " << ft2_ntotal << " entries, of which " 
	  << ft2exp.nInt() << " overlap with GTI\n";
      
      if(!exp_blob_fn.empty())
	{
	  BLOBSerializer blob(exp_blob_fn);
	  if(verbose>=2)
	    std::cout << "Saving exposure BLOB: " 
		      << exp_blob_fn << '\n';
	  ft2exp.partiallySerializeToBlob(blob);
	}
    }

  if(verbose>=2)
    {
      std::cout
	<< "FT2: Total time exposure: " 
	<< lrint(ft2exp.tInt().back().livetime_int_stop)
	<< " seconds\n";
      for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
	std::cout
	  << "FT2: Total exposure (EC=" << iea << "): " 
	  << ft2exp.eInt(iea).back().stop << " cm^2 seconds\n";
    }

  if(!day_exp_lc_fn.empty())
    {
      const std::vector<std::pair<unsigned,double> >& lc(ft2exp.getDayExpLC());
      if(verbose>=2)
	std::cout
	  << "FT2: writing exposure lightcurve to: " << day_exp_lc_fn << '\n';
      std::ofstream str(day_exp_lc_fn.c_str());
      for(unsigned ilc=0;ilc<lc.size();ilc++)
	str << lc[ilc].first << ' ' << lc[ilc].second << '\n';
    }

  if(0)
    {
      std::ofstream str("exposure_dump.dat");
      for(unsigned iint = 0;iint<ft2exp.nInt();iint++)
	{
	  str << std::setprecision(15) << ft2exp.tInt()[iint].start << ' '
	      << std::setprecision(15) << ft2exp.tInt()[iint].stop;
	  for(unsigned iea=0;iea<ft2exp.nEInt();iea++)
	    str << ' ' << ft2exp.eInt(iea)[iint].dexp;
	  str << '\n';
	}
    }

  // --------------------------------------------------------------------------
  // READ IN FT1
  // --------------------------------------------------------------------------

  FT1ROIFilter ft1(ft2exp);
  if(has_roi_cuts && !has_psf_roi_cuts)ft1.setROICuts(ra, dec, roi);
  if(has_psf_roi_cuts)
    ft1.setPSFBasedROICuts(ra, dec, psf_roi_cut_frac, irfs);
  if(has_energy_cuts)ft1.setEnergyCuts(emin, emax);
  FITSVectorDispatcher<FT1> ft1_dispatcher(&ft1);
  FT1::FITSFillOptions fill_opt;
  if(!src_prob_name.empty())fill_opt = FT1::FITSFillOptions(src_prob_name);
  if(verbose>=2)
    std::cout << "Load: FT1 - " << ft1_fn << '\n';
  unsigned ft1_ntotal = 
    ft1_dispatcher.dispatchVector(ft1_fn, FT1::tableName(), fill_opt);
  std::sort(ft1.events.begin(),ft1.events.end());

  if(verbose>=2)
    {
      std:: cout
	<< "FT1: Found " << ft1_ntotal << " events, kept " 
	<< ft1.events.size() << " (";
      for(unsigned ict=0;ict<ft1.nCT();ict++)
	{
	  if(ict!=0)std::cout << ", ";
	  std::cout << "CT=" << ict << ":" << ft1.nEvents(ict);
	}
      std::cout	<< ")\n";
      if(ft1.nCutNoFT2())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutNoFT2() << " no FT2 entry\n";
      if(ft1.nCutZeroExposure())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutZeroExposure() << " zero exposure\n";
      if(ft1.nCutROI())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutROI() << " outside ROI\n";
      if(ft1.nCutEnergy())
	std::cout 
	  << "FT1: Cuts: " << ft1.nCutEnergy() << " outside energy range\n";
    }

  // --------------------------------------------------------------------------
  // COMPUTE PERIODOGRAM
  // --------------------------------------------------------------------------

  if(verbose>=2)
    {
      std::cout << "Computing periodogram...\n";
    }

  try
    {
      HSLike::Mode mode = HSLike::HMM_SRC_ONLY;
      if(!src_prob_name.empty())mode = HSLike::HMM_SRCBKG_PROB;
      Scanner scanner(scan0, scanN, nscan, nharmonic, ft2exp, 
		      verbose, log_scan, mode);
      scanner.scan(ft1.events,known_freq);
    }
  catch(const std::string& x)
    {
      std::cout << "Caught: " << x << '\n';
      throw;
    }
  
  // --------------------------------------------------------------------------
  // CLEAN UP
  // --------------------------------------------------------------------------

  for(std::vector<irfInterface::Irfs *>::iterator iirf=irfs.begin();
      iirf!=irfs.end();iirf++)delete *iirf;
  for(std::vector<PLIntegratedEA*>::iterator ie = eaint.begin();
      ie!=eaint.end(); ie++)delete *ie;
}
