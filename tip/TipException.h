/** \file TipException.h

    \brief Exceptions thrown by table objects.

    \author James Peachey, HEASARC
*/
#ifndef tip_TipException_h
#define tip_TipException_h

#include <exception>
#include <string>

namespace tip {

  /** \class TipException

      \brief Base class for exceptions thrown by table objects.
  */
  class TipException : public std::exception {
    public:
      TipException(const std::string & msg = "Table component exception"): std::exception(), m_msg(msg), m_status(1) {}
      TipException(int status, const std::string & msg = "Table component exception");
      virtual ~TipException() throw() {}
      virtual const char * what() const throw() { return m_msg.c_str(); }
      virtual int code() const { return m_status; }

    protected:
      std::string m_msg;
      int m_status;
  };

}

#endif
