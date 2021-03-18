/** \file FileSummary.cxx
    \brief Data file-wide access.
    \author James Peachey, HEASARC
*/
#include <string>

#include "tip/FileSummary.h"

namespace tip {

  ExtSummary::ExtSummary(const std::string & ext_id): m_ext_id(ext_id) {}

  const std::string & ExtSummary::getExtId() const { return m_ext_id; }
}
