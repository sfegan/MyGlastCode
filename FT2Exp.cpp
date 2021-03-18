//-*-mode:c++; mode:font-lock;-*-

#include "VSAAlgebra.hpp"
#include "FT2Exp.hpp"
#include "Util.hpp"

using namespace VERITAS::VSAAlgebra;

FT2Exp::~FT2Exp()
{
  // nothing to see here
}

void FT2Exp::visitAcceptedInterval(unsigned irow, FT2& ft2)
{
  if(!m_t0set)m_t0=ft2.t_start, m_t0set=true;

  double dt = ft2.livetime;
  double livetime_frac = dt/(ft2.t_stop-ft2.t_start);

  // Replace with Quaternion method since we need phi also
  // const double ct = 
  //   std::cos(sphere_dist(m_ra, m_dec, ft2.scz_ra, ft2.scz_dec));
  
  RotationVec3D sc2eq = ft2.vecScRot();
  Vec3D src_vec = Vec3D::makeLatLon(m_dec,m_ra);
  src_vec.rotate(-sc2eq);
  const double ct = src_vec.cosTheta();
  const double phi = src_vec.phi();

  m_tacc.add(dt);
  m_tint.push_back(TimeInterval(ft2.t_start-m_t0,ft2.t_stop-m_t0,dt,
				m_tacc.sum(),ct,phi));

  unsigned bin = lrint(ft2.t_start/m_lcperiod-0.5);
  if(bin != m_lcbin)
    {
      if(m_lcbin != 0)m_lcexp.push_back(std::make_pair<double,double>
					 ((double(m_lcbin)+0.5)*m_lcperiod,
					  m_lcacc.sum()));
      m_lcacc.reset();
      m_lcbin = bin;
    }

  for(unsigned iea=0;iea<m_ea.size();iea++)
    if(m_ea[iea])
      {
	double aeff = m_ea[iea]->value(ct, phi, livetime_frac)*10000;

	ExposureInterval i;
	i.exp_int_start               = m_eacc[iea].sum();
	i.aeff                        = aeff;
	double exposure = aeff*dt;
	i.dexp                        = exposure;
	m_eacc[iea].add(exposure);
	i.exp_int_stop                = m_eacc[iea].sum();
	m_eint[iea].push_back(i);
	m_lcacc.add(exposure);
      }
}

void FT2Exp::leaveFileVector()
{
  if(m_lcbin != 0)m_lcexp.push_back(std::make_pair<double,double>
				    ((double(m_lcbin)+0.5)*m_lcperiod,
				     m_lcacc.sum()));
}

unsigned FT2Exp::findTInterval(double t) const
{
  return _findInterval(t,m_tint);
}

unsigned FT2Exp::findEInterval(double t, unsigned ie) const
{
  return _findInterval(t,m_eint[ie]);
}

void FT2Exp::mergeExposures()
{
  const unsigned nint = m_tint.size();
  m_eint[0].resize(nint);
  for(unsigned ie=1;ie<m_eint.size();ie++)
    if(m_eint[ie].empty())
      for(unsigned iint=0;iint<nint;iint++)
	{
	  m_eint[0][iint].exp_int_start += m_eint[ie][iint].exp_int_start;
	  m_eint[0][iint].exp_int_stop  += m_eint[ie][iint].exp_int_stop;
	  m_eint[0][iint].dexp          += m_eint[ie][iint].dexp;
	  m_eint[0][iint].aeff          += m_eint[ie][iint].aeff;
	}
  m_eint.resize(1);
  m_merged_exposures=true;
}

bool FT2Exp::partiallySerializeToBlob(BLOBSerializer& s) const
{
  return FT2ROI::serializeToBlob(s)
    && s.serialize(m_t0set)
    && s.serialize(m_t0)
    && s.serialize(m_tint)
    && s.serialize(m_eint)
    && s.serialize(m_merged_exposures)
    && s.serialize(m_lcperiod)
    && s.serialize(m_lcexp);
}

bool FT2Exp::partiallyUnserializeFromBlob(BLOBUnserializer& s)
{
  return FT2ROI::unserializeFromBlob(s)
    && s.unserialize(m_t0set)
    && s.unserialize(m_t0)
    && s.unserialize(m_tint)
    && s.unserialize(m_eint)
    && s.unserialize(m_merged_exposures)
    && s.unserialize(m_lcperiod)
    && s.unserialize(m_lcexp);
}
