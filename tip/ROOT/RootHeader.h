/** \file RootHeader.h

    \brief High level encapsulation of a FITS-like header.

    \author James Peachey, HEASARC
*/
#ifndef tip_RootHeader_h
#define tip_RootHeader_h

#include <string>

#include "tip/Header.h"

namespace tip {

  /** \class RootHeader

      \brief High level encapsulation of a FITS-like header.
  */
  class RootHeader : public Header {
    public:
      /** \brief Return the name of the particular column implementation (subclass identifier).
      */
      virtual const std::string implementation() const { return "Root"; }

  };

}

#endif
