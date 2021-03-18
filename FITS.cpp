#include "FITS.hpp"

// ----------------------------------------------------------------------------
// 
// FITS header
//
// ----------------------------------------------------------------------------

void FITSHeader::loadFromFITS(const std::string& filename, 
			      const std::string& tablename)
{
  const tip::Table* fits =
    tip::IFileSvc::instance().readTable(filename, tablename);

  const tip::Header & header(fits->getHeader());

  loadFromHeader(header);
}

void FITSHeader::loadFromHeader(const tip::Header & header)
{
  tip::Header::ConstIterator ih = header.begin();
  while(ih != header.end())
    {
      m_kv[ih->getName()].push_back(ih->getValue());
      ih++;
    }
}
