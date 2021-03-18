/** \file Extension.h

    \brief High level encapsulation of a FITS-like file extension, which may or may not actually BE in FITS format.

    \author James Peachey, HEASARC
*/
#ifndef tip_Extension_h
#define tip_Extension_h

namespace tip {

  class Header;

  /** \class Extension

      \brief High level encapsulation of a FITS-like file extension.
  */
  class Extension {
    public:
      /** \brief Destruct an extension object.
      */
      virtual ~Extension() {}

      /** \brief Retrieve Header object, which is a container of FITS-like keywords, non-const version.
      */
      virtual Header & getHeader() = 0;

      /** \brief Retrieve Header object, which is a container of FITS-like keywords, const version.
      */
      virtual const Header & getHeader() const = 0;

      /** \brief Returns true if the extension is a image, false otherwise.
      */
      virtual bool isImage() const = 0;

      /** \brief Returns true if the extension is a table, false otherwise.
      */
      virtual bool isTable() const = 0;

      /** \brief Return name of this extension.
      */
      virtual const std::string & getName() const = 0;

      /** \brief Set name of this extension.
      */
      virtual void setName(const std::string & name) = 0;
  };

}

#endif
