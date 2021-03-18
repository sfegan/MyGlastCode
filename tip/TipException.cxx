/** \file TipException.cxx

    \brief Exceptions thrown by table objects.

    \author James Peachey, HEASARC
*/
#include <sstream>
#include <string>

#include "fitsio.h"
#include "tip/TipException.h"

namespace tip {

  TipException::TipException(int status, const std::string & msg): std::exception(), m_msg(), m_status(status) {
    static const std::string fitsio_unknown_message("unknown error status");
    char fitsio_error_msg[31] = "";
    fits_get_errstatus(m_status, fitsio_error_msg);
    if (fitsio_error_msg == fitsio_unknown_message) {
      m_msg = msg;
    } else {
      std::ostringstream os;
      os << msg << " (CFITSIO ERROR " << m_status << ": " << fitsio_error_msg << ")";
      m_msg = os.str();
    }
  }
}

