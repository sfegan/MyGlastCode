/** \file Header.h

    \brief High level encapsulation of a FITS-like header.

    \author James Peachey, HEASARC
*/
#ifndef tip_Header_h
#define tip_Header_h

#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "tip/Iterator.h"
#include "tip/KeyRecord.h"
#include "tip/TipException.h"

namespace tip {

  class Header;

  /** \class Keyword

      \brief Encapsulation of a single keyword. The keyword may be get/set in any data type.

             See also the KeyRecord class. The difference between them is that Keyword obejcts
             are bound to a specific data header, whereas KeyRecord objects simply contain a
             copy of the entire keyword record stored as a single string.
  */
  class Keyword {
    public:
      static Keyword & emptyKeyword();

      /** \brief Construct a Keyword object associated with a particular Header.
          \param header_data Pointer to the referent IExtensionData object.
          \param name The name of this Keyword.
      */
      Keyword(Header * header_data, const std::string & name): m_record(), m_header_data(header_data),
        m_name(name) {
        if (0 == m_header_data) throw TipException("Keyword::Keyword(Header *, const std::string &): "
          "Cannot create Keyword with a NULL IExtensionData pointer");
      }

      /** \brief Get the current value of this Keyword.
          The type of the converted value is given by the template parameter.
          \param value The current value.
      */
      template <typename T>
      void get(T & value) const;

      /** \brief Set the value of this Keyword.
          \param value The value to set.
      */
      template <typename T>
      void set(const T & value);

      void getRecord(KeyRecord & rec) const;

      void setRecord(const KeyRecord & rec);

      /// \brief Get comment associated with this keyword.
      std::string getComment() const;

      /// \brief Set comment associated with this keyword.
      void setComment(const std::string & comment);

      /// \brief Get unit associated with this keyword.
      std::string getUnit() const;

      /** \brief Set unit associated with this keyword.
          \param unit The new unit of the keyword.
      */
      void setUnit(const std::string & unit);

    private:
      Keyword(): m_record(), m_header_data(0), m_name() {}

      KeyRecord m_record;
      Header * m_header_data;
      std::string m_name;
  };

  /** \class Header

      \brief High level encapsulation of a FITS-like header.
  */
  class Header {
    public:
      /** \brief For convenience typedef the underlying keyword container.
      */
      typedef std::map<std::string, Keyword> KeywordCont_t;

      /** \brief Container of keywords in order.
      */
      typedef std::vector<KeyRecord> KeySeq_t;
      typedef KeySeq_t::iterator Iterator;
      typedef KeySeq_t::const_iterator ConstIterator;

      /** \brief Adaptor for keywords in a form separate from the Header structure.
      */
      typedef std::pair<std::string, std::string> KeyValPair_t;

      /** \brief Container of key-value pairs.
      */
      typedef std::vector<KeyValPair_t> KeyValCont_t;

      virtual ~Header() {}

      /** \brief Random read/write keyword access.
          \param name The name of the keyword.
      */
      Keyword & operator [](const std::string & name) { return find_or_make(name); }

      /** \brief Random read only keyword access.
          \param name The name of the keyword.
      */
      const Keyword & operator [](const std::string & name) const { return find_or_make(name); }

      /** \brief Obtain a list of keywords in a container of key-value pairs.
          \param keys The list of keys to obtain. MUST BE NULL TERMINATED!
          \param key_vals The output container of key-value pairs
      */
      void get(const char ** keys, KeyValCont_t & key_vals) const;

      /** \brief Update header's keywords using the container of key-value pairs.
          New keywords are not added; only keywords which already exist will be updated.
          \param key_vals The input container of key-value pairs
      */
      void update(const KeyValCont_t & key_vals);

      /** \brief Return a time in the standard format for use as a keyword.
          \param time The time to format.
      */
      std::string formatTime(const time_t & time) const;

      virtual KeySeq_t::size_type getNumKeywords() const { unsupported("getNumKeywords()"); return 0; }

      virtual Iterator begin() { unsupported("begin()"); return KeySeq_t().begin(); }

      virtual Iterator end() { unsupported("end()"); return KeySeq_t().end(); }

      virtual ConstIterator begin() const { unsupported("begin() const"); return KeySeq_t().begin(); }

      virtual ConstIterator end() const { unsupported("end() const"); return KeySeq_t().end(); }

      /** \brief Return an iterator pointing to the first keyword with the given name. If no keyword
          with this name was found, returns end().
          \param key_name The name of the keyword being sought.
      */
      virtual Iterator find(const std::string &) { unsupported("find(const std::string &)"); return KeySeq_t().end(); }

      /** \brief Return a const iterator pointing to the first keyword with the given name. If no keyword
          with this name was found, returns end().
          \param key_name The name of the keyword being sought.
      */
      virtual ConstIterator find(const std::string &) const
        { unsupported("find(const std::string &) const"); return KeySeq_t().end(); }

      /** \brief Insert a keyword record before the given iterator position.
          with this name was found, returns end().
          \param itor The position before which the new record will be placed.
          \param record The record being inserted.
      */
      virtual Iterator insert(Iterator itor, const KeyRecord &) { unsupported("insert(Iterator, const KeyRecord &)"); return itor; }

      /** \brief Append a keyword record to the container of keywords in this header.
          This has the same effect as calling insert(end(), record).
          \param record The record being appeneded.
      */
      virtual Iterator append(const KeyRecord &) { unsupported("append(const KeyRecord &)"); return KeySeq_t().end(); }

      /** \brief Erase the keyword pointed to by the iterator. Other keywords with the same name will
          not be erased.
          \param itor Iterator pointing to the keyword being erased.
      */
      virtual Iterator erase(Iterator) { unsupported("erase(Iterator)"); return KeySeq_t().end(); }

      /** \brief Erase all keywords whose name is the same as the given name.
          \param key_name The name of the keyword(s) being erased.
      */
      virtual void erase(const std::string &) { unsupported("erase(const std::string &)"); }

      /** \brief Get a keyword from this header data object.
          \param name The name of the keyword to get from the header data object.
          \param value The output value of the keyword, converted to the given type.
      */
      virtual void getKeyword(const std::string &, bool &) const
        { unsupported("getKeyword(const std::string &, bool &) const"); }
      virtual void getKeyword(const std::string &, double &) const
        { unsupported("getKeyword(const std::string &, double &) const"); }
      virtual void getKeyword(const std::string &, float &) const
        { unsupported("getKeyword(const std::string &, float &) const"); }
      virtual void getKeyword(const std::string &, char &) const
        { unsupported("getKeyword(const std::string &, char &) const"); }
      virtual void getKeyword(const std::string &, signed char &) const
        { unsupported("getKeyword(const std::string &, signed char &) const"); }
      virtual void getKeyword(const std::string &, signed short &) const
        { unsupported("getKeyword(const std::string &, signed short &) const"); }
      virtual void getKeyword(const std::string &, signed int &) const
        { unsupported("getKeyword(const std::string &, signed int &) const"); }
      virtual void getKeyword(const std::string &, signed long &) const
        { unsupported("getKeyword(const std::string &, signed long &) const"); }
      virtual void getKeyword(const std::string &, unsigned char &) const
        { unsupported("getKeyword(const std::string &, unsigned char &) const"); }
      virtual void getKeyword(const std::string &, unsigned short &) const
        { unsupported("getKeyword(const std::string &, unsigned short &) const"); }
      virtual void getKeyword(const std::string &, unsigned int &) const
        { unsupported("getKeyword(const std::string &, unsigned int &) const"); }
      virtual void getKeyword(const std::string &, unsigned long &) const
        { unsupported("getKeyword(const std::string &, unsigned long &) const"); }
      virtual void getKeyword(const std::string &, std::string &) const
        { unsupported("getKeyword(const std::string &, std::string &) const"); }

      virtual void getKeyRecord(const std::string &, std::string &) const
        { unsupported("getKeyRecord(const std::string &, std::string &) const"); }

      /** \brief Set a keyword in this header data object.
          \param name The name of the keyword to set in the header data object.
          \param value The input value of the keyword.
      */
      virtual void setKeyword(const std::string &, const bool &)
        { unsupported("setKeyword(const std::string &, const bool &)"); }
      virtual void setKeyword(const std::string &, const double &)
        { unsupported("setKeyword(const std::string &, const double &)"); }
      virtual void setKeyword(const std::string &, const float &)
        { unsupported("setKeyword(const std::string &, const float &)"); }
      virtual void setKeyword(const std::string &, const char &)
        { unsupported("setKeyword(const std::string &, const char &)"); }
      virtual void setKeyword(const std::string &, const signed char &)
        { unsupported("setKeyword(const std::string &, const signed char &)"); }
      virtual void setKeyword(const std::string &, const signed short &)
        { unsupported("setKeyword(const std::string &, const signed short &)"); }
      virtual void setKeyword(const std::string &, const signed int &)
        { unsupported("setKeyword(const std::string &, const signed int &)"); }
      virtual void setKeyword(const std::string &, const signed long &)
        { unsupported("setKeyword(const std::string &, const signed long &)"); }
      virtual void setKeyword(const std::string &, const unsigned char &)
        { unsupported("setKeyword(const std::string &, const unsigned char &)"); }
      virtual void setKeyword(const std::string &, const unsigned short &)
        { unsupported("setKeyword(const std::string &, const unsigned short &)"); }
      virtual void setKeyword(const std::string &, const unsigned int &)
        { unsupported("setKeyword(const std::string &, const unsigned int &)"); }
      virtual void setKeyword(const std::string &, const unsigned long &)
        { unsupported("setKeyword(const std::string &, const unsigned long &)"); }
      virtual void setKeyword(const std::string &, const std::string &)
        { unsupported("setKeyword(const std::string &, const std::string &)"); }
      virtual void setKeyword(const std::string &, const char * const &)
        { unsupported("setKeyword(const std::string &, const char * const &)"); }

      virtual void setKeyRecord(const std::string &, const std::string &)
        { unsupported("setKeyRecord(const std::string &, const std::string &)"); }

      /** \brief Return the name of the particular column implementation (subclass identifier).
      */
      virtual const std::string implementation() const = 0;

      /** \brief Get comment associated with the given keyword.
          \param name The name of the keyword.
      */
      virtual std::string getKeyComment(const std::string &) const
        { unsupported("getKeyComment(const std::string &)"); return ""; } 

      /** \brief Set comment associated with the given keyword.
          \param name The name of the keyword.
          \param comment The new comment of the keyword.
      */
      virtual void setKeyComment(const std::string &, const std::string &)
        { unsupported("setKeyComment(const std::string &, const std::string &)"); } 

      /** \brief Get unit associated with the given keyword.
          \param name The name of the keyword.
      */
      virtual std::string getKeyUnit(const std::string &) const
        { unsupported("getKeyUnit(const std::string &)"); return ""; } 

      /** \brief Set unit associated with the given keyword.
          \param name The name of the keyword.
          \param unit The new unit of the keyword.
      */
      virtual void setKeyUnit(const std::string &, const std::string &)
        { unsupported("setKeyUnit(const std::string &, const std::string &)"); } 

      /** \brief Add a descriptive comment message to the header.
          \param comment The comment string to add.
      */
      virtual void addComment(const std::string & /* comment */)
        { unsupported("addComment(const std::string & comment)"); }

      /** \brief Add a descriptive history message to the header.
          \param history The history string to add.
      */
      virtual void addHistory(const std::string & /* history */)
        { unsupported("addHistory(const std::string & history)"); }

    protected:
      /** \brief Internal utility to add keywords when they are looked up.
      */
      Keyword & find_or_make(const std::string & name) const;

    private:
      void unsupported(const std::string & method) const {
        throw TipException(std::string("Header method ") + method + " is not supported for the " + implementation() +
          " implementation");
      }

      KeywordCont_t m_keywords;
  };

  inline Keyword & Header::find_or_make(const std::string & name) const {
    // Because inquiries for keywords for a constant object can still add them to the Header, need to
    // get rid of const.
    Header & header = const_cast<Header &>(*this);

    // Look for keyword.
    KeywordCont_t::iterator itor = header.m_keywords.find(name);

    // If not found, create a new one.
    if (header.m_keywords.end() == itor)
      itor = header.m_keywords.insert(itor, std::make_pair(name, Keyword(&header, name)));

    // Return the keyword.
    return itor->second;
  }

}

namespace tip {

  template <typename T>
  inline void Keyword::get(T & value) const { m_header_data->getKeyword(m_name, value); }

  template <typename T>
  inline void Keyword::set(const T & value) { m_header_data->setKeyword(m_name, value); }

  inline void Keyword::getRecord(KeyRecord & rec) const {
    std::string rec_str;
    m_header_data->getKeyRecord(m_name, rec_str);
    rec.set(rec_str);
  }

  inline void Keyword::setRecord(const KeyRecord & rec) {
    m_header_data->setKeyRecord(m_name, rec.get());
  }

  inline std::string Keyword::getComment() const {
    return m_header_data->getKeyComment(m_name);
  }

  inline void Keyword::setComment(const std::string & comment) {
    m_header_data->setKeyComment(m_name, comment);
  }

  inline std::string Keyword::getUnit() const {
    return m_header_data->getKeyUnit(m_name);
  }

  inline void Keyword::setUnit(const std::string & unit) {
    m_header_data->setKeyUnit(m_name, unit);
  }

}

#endif
