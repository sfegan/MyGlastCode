//-*-mode:c++; mode:font-lock;-*-

#include "Util.hpp"
#include "FT1ROIFilter.hpp"

double FT1ROIFilter::
getPSFRegion(unsigned iirf, double energy, double theta, double phi)
{
  irfInterface::IPsf* psf = m_irfs[iirf]->psf();
  double dl = 0;
  double fl = 0;
  double dr = 180;
  double fr = 1.0;
  while((fr-fl)>(1-m_psf_frac)*0.001)
    {
      double dm = 0.5*(dl+dr);
      double fm = psf->angularIntegral(energy, theta, phi, dm, 0);
#if 0
      std::cout << m_psf_frac << ' '
	<< dl << ' ' << dr << ' ' << dm << ' '
	<< fl << ' ' << fr << ' ' << fm << '\n';
#endif
      if(fm > m_psf_frac)dr=dm,fr=fm;
      else dl=dm,fl=fm;
    }
  //  std::cout << "E: " << energy << ' '<< 0.5*(dl+dr)/180.0*M_PI << '\n';
  return 0.5*(dl+dr)/180.0*M_PI;
}

void FT1ROIFilter::visitElement(unsigned irow, FT1& ft1)
{
  std::pair<unsigned,unsigned> eventkey(ft1.run_id,ft1.event_id);
  std::set<std::pair<unsigned,unsigned> >::const_iterator ikey = 
    m_events_seen.find(eventkey);
  if(ikey!=m_events_seen.end())
    { m_ncut_duplicate++; return; }
  m_events_seen.insert(eventkey);
  if((m_energy_cuts)&&(ft1.energy<m_emin || ft1.energy>m_emax))
    { m_ncut_energy++; return; }
  if((m_roi_cuts)&&(sphere_dist(m_ra,m_dec,ft1.ra,ft1.dec)>m_roi))
    { m_ncut_roi++; return; }
  double evtime = ft1.time-m_ft2exp.t0();
  unsigned iint = m_ft2exp.findTInterval(evtime);
  if(iint == m_ft2exp.nInt())
    { m_ncut_ft2exp++; return; }

  unsigned itype = unsigned(ft1.conversion_type);
  if(m_ft2exp.hasMergedExposures())itype=0;

  if((itype >= m_ft2exp.nEInt())||(m_ft2exp.eInt(itype).empty())
     ||(m_ft2exp.eInt(itype)[iint].stop==m_ft2exp.eInt(itype)[iint].start))
    { m_ncut_zeroexp++; return; }

  if(m_psf_roi_cuts)
    { 
      const double theta = std::acos(m_ft2exp.tInt()[iint].costheta)/180*M_PI;
      const double phi = m_ft2exp.tInt()[iint].phi/M_PI*180.0;
      const double d = 
	getPSFRegion(unsigned(ft1.conversion_type),ft1.energy,theta,phi);
      if(sphere_dist(m_ra,m_dec,ft1.ra,ft1.dec)>d)
	{ m_ncut_roi++; return; }
    }
  
  if(m_nevent.size()<=itype)m_nevent.resize(itype+1);
  m_nevent[itype]++;

  events.push_back(FT1AEff(ft1, m_ft2exp.eInt(itype)[iint].aeff, iint));
}

void SimpleEventTimesSimulation::
generateEventTimes(std::vector<FT1AEff>& events)
{
  events.clear();

  double nest = 0;
  for(unsigned ie=0;ie!=m_ft2exp.nEInt();ie++)
    if(!m_ft2exp.eInt(ie).empty() && ie<m_rates.size() && m_rates[ie]>0)
      nest += m_ft2exp.eInt(ie).back().stop * m_rates[ie];
  events.reserve(lrint(nest*2.0*std::sqrt(nest)));
  
  unsigned ievent = 1;
  for(unsigned ie=0;ie!=m_ft2exp.nEInt();ie++)
    if(!m_ft2exp.eInt(ie).empty() && ie<m_rates.size() && m_rates[ie]>0)
      {
	const unsigned nint = m_ft2exp.eInt(ie).size();
	const double eend = m_ft2exp.eInt(ie).back().stop;
#if 1
	double e = m_rng.Exponential()/m_rates[ie];
	while(e<eend)
	  {
#else
	const unsigned nev = eend*m_rates[ie];//m_rng.Poisson(eend*m_rates[ie]);
	for(unsigned iev=0;iev<nev;iev++)
	  {
	    const double e = m_rng.Uniform()*eend;
#endif
	    unsigned iint = m_ft2exp.findEInterval(e,ie);
	    assert(iint<nint);
	    const FT2Exp::ExposureInterval& E(m_ft2exp.eInt(ie)[iint]);
	    const FT2Exp::TimeInterval& T(m_ft2exp.tInt()[iint]);
	    double x = (e-E.start)/(E.stop-E.start);
	    FT1AEff ft1;
	    ft1.aeff = E.aeff;
	    ft1.ift2exp = iint;
	    ft1.conversion_type = ie;
	    ft1.time = (1.0-x)*T.start + x*T.stop + m_ft2exp.t0();
	    ft1.event_id = ievent++;
	    ft1.run_id = lrint(ft1.time/3600-0.5)*3600;
#if 0
	    std::cout << ft1.time << ' ' << iint << ' ' << x << ' '
		      << e << ' ' << E.start << ' ' << E.start << '\n';
#endif
	    events.push_back(ft1);
#if 1
	    e += m_rng.Exponential()/m_rates[ie];
#endif
	  }
      }
  std::sort(events.begin(),events.end());
}
