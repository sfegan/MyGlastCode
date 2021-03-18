//-*-mode:c++; mode:font-lock;-*-

#include "FT2ROI.hpp"
#include "Util.hpp"

FT2ROI::~FT2ROI()
{
  // nothing to see here
}

void FT2ROI::visitElement(unsigned irow, FT2& ft2)
{
  if(m_rnlm && ft2.lat_mode == 4)
    {
      visitRejectedInterval(irow,ft2);
      return;
    }

  double livetime_frac = ft2.livetime/(ft2.t_stop-ft2.t_start);

  GTIRange::const_iterator igti;  
  igti = m_ft2gti.findGTI(ft2.t_start);
  if(igti!=m_ft2gti.end())
    {
      if(igti->t_stop<ft2.t_stop)
	{
	  ft2.t_start=igti->t_stop;
	  ft2.livetime = livetime_frac*(ft2.t_stop-ft2.t_start);
	}
      else
	{
	  visitRejectedInterval(irow,ft2);
	  return;
	}
    }

  igti = m_ft2gti.findGTI(ft2.t_stop);
  if(igti!=m_ft2gti.end())
    {
      if(igti->t_start>ft2.t_start)
	{
	  ft2.t_stop=igti->t_start;
	  ft2.livetime = livetime_frac*(ft2.t_stop-ft2.t_start);
	}
      else
	{
	  visitRejectedInterval(irow,ft2);
	  return;
	}
    }  

  igti = m_gti.findGTI(ft2.t_start);
  if(igti==m_gti.end())
    {
      igti = m_gti.findGTI(ft2.t_stop);
      if(igti==m_gti.end())
	{
	  visitRejectedInterval(irow,ft2);
	  return;
	}
      ft2.t_start = igti->t_start;
    }
  else if(ft2.t_stop>igti->t_stop)
    {
      ft2.t_stop=igti->t_stop;
    }

  if(ft2.t_stop<=ft2.t_start)
    {
      visitRejectedInterval(irow,ft2);
      return;
    }

  if(m_zmax < M_PI)
    {
      const double zroi = sphere_dist(m_ra, m_dec, ft2.zn_ra, ft2.zn_dec);
      if(zroi>m_zmax)
	{
	  visitRejectedInterval(irow,ft2);
	  return;
	}
    }

  if(m_tmin > 0 && m_tmax < M_PI)
    {
      const double troi = sphere_dist(m_ra, m_dec, ft2.scz_ra, ft2.scz_dec);
      if(troi<m_tmin || troi>m_tmax)
	{
	  visitRejectedInterval(irow,ft2);
	  return;
	}
    }

  m_ft2gti.addGTI(ft2);

#if 0
  std::cout << ft2.in_saa << ' ' <<  ft2.lat_mode << ' '
	    << ft2.lat_config << ' ' << ft2.data_qual << '\n';
#endif

  visitAcceptedInterval(irow,ft2);
}

void FT2ROI::visitAcceptedInterval(unsigned irow, FT2& ft2)
{
  // nothing to see here
}

void FT2ROI::visitRejectedInterval(unsigned irow, FT2& ft2)
{
  // nothing to see here
}

bool FT2ROI::serializeToBlob(BLOBSerializer& s) const
{
  return s.serialize(m_gti)
    && s.serialize(m_ft2gti)
    && s.serialize(m_ra)
    && s.serialize(m_dec)
    && s.serialize(m_zmax)
    && s.serialize(m_tmin)
    && s.serialize(m_tmax)
    && s.serialize(m_rnlm);
}

bool FT2ROI::unserializeFromBlob(BLOBUnserializer& s)
{
  return s.unserialize(m_gti)
    && s.unserialize(m_ft2gti)
    && s.unserialize(m_ra)
    && s.unserialize(m_dec)
    && s.unserialize(m_zmax)
    && s.unserialize(m_tmin)
    && s.unserialize(m_tmax)
    && s.unserialize(m_rnlm);
}
