//-*-mode:c++; mode:font-lock;-*-

#include <algorithm>

#include "Accumulator.hpp"
#include "GTI.hpp"

bool GTI::fillFromFITS(tip::Table::ConstRecord & sc_datum, 
		       const FITSHeader& header,
		       const FITSFillOptions& opt)
{
  sc_datum["start"].get(t_start);
  sc_datum["stop"].get(t_stop);

  return true;
}

unsigned GTI::loadFromFITS(std::vector<GTI>& GTI_vec,
			   const std::string& filename, 
			   const std::string& tablename,
			   const FITSFillOptions& opt)
{
  return fillVectorFromFITS<GTI>(GTI_vec, filename, tablename, opt);  
}

double GTIRange::totalTime() const
{
  Accumulator a;
  for(std::vector<GTI>::const_iterator igti=m_gti.begin();
      igti!=m_gti.end();igti++)a.add(igti->t_stop-igti->t_start);
  return a.sum();
}

void GTIRange::loadGTIsFromFITS(const std::string& filename,   
				const std::string& tablename,
				const GTI::FITSFillOptions& opt)
{
  GTI::loadFromFITS(m_gti, filename, tablename, opt);
  normalize();
}

void GTIRange::normalize()
{
  std::sort(m_gti.begin(), m_gti.end());
  unsigned jgti=0;
  for(unsigned igti=1;igti<m_gti.size();igti++)
    {
      if(m_gti[igti].t_start <= m_gti[jgti].t_stop)
	{
	  // Overlap: extend the previous interval if necessary
	  if(m_gti[igti].t_stop > m_gti[jgti].t_stop)
	    m_gti[jgti].t_stop = m_gti[igti].t_stop;
	}
      else
	{
	  jgti++;
	  if(igti != jgti)m_gti[jgti]=m_gti[igti];
	}
    }
  jgti++;
  m_gti.resize(jgti);
  m_last = begin();
}

bool GTIRange::serializeToBlob(BLOBSerializer& s) const
{ 
  return s.serialize(m_gti); 
}

bool GTIRange::unserializeFromBlob(BLOBUnserializer& s)
{ 
  bool ok = s.unserialize(m_gti); 
  m_last=m_gti.begin(); 
  return ok;
}
