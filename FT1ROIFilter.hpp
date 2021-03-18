//-*-mode:c++; mode:font-lock;-*-

#ifndef FT1ROIFILTER_HPP
#define FT1ROIFILTER_HPP

#include <RandomNumbers.hpp>

#include <irfInterface/IEfficiencyFactor.h>
#include <irfInterface/IrfsFactory.h>
#include <irfLoader/Loader.h>

#include "FT1.hpp"
#include "FT2Exp.hpp"

class FT1AEff: public FT1
{
public:
  FT1AEff(): FT1(), aeff(), ift2exp() { }
  FT1AEff(const FT1& ft1, double _aeff, unsigned iexp): 
    FT1(ft1), aeff(_aeff), ift2exp(iexp) { }
  double aeff;
  unsigned ift2exp;
};

inline bool operator<(const FT1AEff& a, const FT1AEff& b) 
{
  return a.time<b.time;
}

class FT1ROIFilter: public FITSVectorVisitor<FT1>
{
public:
  FT1ROIFilter(const FT2Exp& ft2exp): 
    FITSVectorVisitor<FT1>(),
    events(), 
    m_ft2exp(ft2exp), m_ncut_ft2exp(), m_ncut_zeroexp(), 
    m_roi_cuts(), m_psf_roi_cuts(),
    m_ra(), m_dec(), m_roi(), m_psf_frac(), m_irfs(), m_ncut_roi(),
    m_energy_cuts(), m_emin(), m_emax(), m_ncut_energy(),
    m_nevent(), m_events_seen(), m_ncut_duplicate() { }

  void setROICuts(double ra, double dec, double roi)
  { m_roi_cuts=true; m_ra=ra; m_dec=dec; m_roi=roi; }
  void setPSFBasedROICuts(double ra, double dec, double psf_frac,
			  std::vector<irfInterface::Irfs *> irfs)
  { m_psf_roi_cuts=true; m_ra=ra; m_dec=dec; m_psf_frac=psf_frac; m_irfs=irfs;}
  void setEnergyCuts(double emin, double emax) 
  { m_energy_cuts=true; m_emin=emin; m_emax=emax; }

  virtual void visitElement(unsigned irow, FT1& ft1);

  unsigned nCutNoFT2() const { return m_ncut_ft2exp; }
  unsigned nCutZeroExposure() const { return m_ncut_zeroexp; }
  unsigned nCutROI() const { return m_ncut_roi; }
  unsigned nCutEnergy() const { return m_ncut_energy; }
  unsigned nCutDuplicate() const { return m_ncut_duplicate; }

  unsigned nEvents(unsigned itype) const
  { if(itype<m_nevent.size())return m_nevent[itype]; else return 0; }
  unsigned nCT() const { return m_nevent.size(); }

  std::vector<FT1AEff> events;
private:  
  double getPSFRegion(unsigned iirf, double energy, double theta, double phi);

  const FT2Exp& m_ft2exp;
  unsigned m_ncut_ft2exp;
  unsigned m_ncut_zeroexp;
  bool m_roi_cuts;
  bool m_psf_roi_cuts;
  double m_ra;
  double m_dec;
  double m_roi;
  double m_psf_frac;
  std::vector<irfInterface::Irfs *> m_irfs;
  unsigned m_ncut_roi;
  bool m_energy_cuts;
  double m_emin;
  double m_emax;
  unsigned m_ncut_energy;
  std::vector<unsigned> m_nevent;
  std::set<std::pair<unsigned,unsigned> > m_events_seen;
  unsigned m_ncut_duplicate;
};

class SimpleEventTimesSimulation
{
public:
  SimpleEventTimesSimulation(RandomNumbers& rng, const FT2Exp& ft2exp, 
			     const std::vector<double>& rates): 
    m_rng(rng), m_ft2exp(ft2exp), m_rates(rates) { }
  void generateEventTimes(std::vector<FT1AEff>& events);
private:
  RandomNumbers&       m_rng;
  const FT2Exp&        m_ft2exp;
  std::vector<double>  m_rates;
};

#endif // #ifndef FT1ROIFILTER_HPP
