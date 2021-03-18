/** \file IColumn.h
    \brief Generic abstract interface to columnar data.
    \author James Peachey, HEASARC
*/
#ifndef tip_IColumn_h
#define tip_IColumn_h

#include <string>
#include <vector>

#include "tip/Header.h"
#include "tip/TipException.h"
#include "tip/tip_types.h"

namespace tip {

  /** \class IColumn
      \brief Generic abstract interface to columnar data.
  */
  class IColumn {
    public:
      /** \brief Construct a column with the given identifier.
      */
      IColumn(const std::string & id = std::string()): m_units(), m_id(id) {}

      virtual ~IColumn() throw() {}

      /** \brief Get a value from this column. An implementation which just throws an exception is
          provided in the base class to simplify derived classes which need only a small number of
          the methods to function.
      */
      virtual void get(Index_t, bool &) const { unsupported("get(Index_t, bool &)"); }
      virtual void get(Index_t, double &) const { unsupported("get(Index_t, double &)"); }
      virtual void get(Index_t, float &) const { unsupported("get(Index_t, float &)"); }
      virtual void get(Index_t, char &) const { unsupported("get(Index_t, char &)"); }
      virtual void get(Index_t, signed char &) const { unsupported("get(Index_t, signed char &)"); }
      virtual void get(Index_t, signed short &) const { unsupported("get(Index_t, signed short &)"); }
      virtual void get(Index_t, signed int &) const { unsupported("get(Index_t, signed int &)"); }
      virtual void get(Index_t, signed long &) const { unsupported("get(Index_t, signed long &)"); }
      virtual void get(Index_t, unsigned char &) const { unsupported("get(Index_t, unsigned char &)"); }
      virtual void get(Index_t, unsigned short &) const { unsupported("get(Index_t, unsigned short &)"); }
      virtual void get(Index_t, unsigned int &) const { unsupported("get(Index_t, unsigned int &)"); }
      virtual void get(Index_t, unsigned long &) const { unsupported("get(Index_t, unsigned long &)"); }
      virtual void get(Index_t, std::string &) const { unsupported("get(Index_t, std::string &)"); }

      virtual void get(Index_t, std::vector<bool> &) const { unsupported("get(Index_t, std::vector<bool> &)"); }
      virtual void get(Index_t, std::vector<double> &) const { unsupported("get(Index_t, std::vector<double> &)"); }
      virtual void get(Index_t, std::vector<float> &) const { unsupported("get(Index_t, std::vector<float> &)"); }
      virtual void get(Index_t, std::vector<char> &) const { unsupported("get(Index_t, std::vector<char> &)"); }
      virtual void get(Index_t, std::vector<signed char> &) const { unsupported("get(Index_t, std::vector<signed char> &)"); }
      virtual void get(Index_t, std::vector<signed short> &) const { unsupported("get(Index_t, std::vector<signed short> &)"); }
      virtual void get(Index_t, std::vector<signed int> &) const { unsupported("get(Index_t, std::vector<signed int> &)"); }
      virtual void get(Index_t, std::vector<signed long> &) const { unsupported("get(Index_t, std::vector<signed long> &)"); }
      virtual void get(Index_t, std::vector<unsigned char> &) const { unsupported("get(Index_t, std::vector<unsigned char> &)"); }
      virtual void get(Index_t, std::vector<unsigned short> &) const { unsupported("get(Index_t, std::vector<unsigned short> &)");}
      virtual void get(Index_t, std::vector<unsigned int> &) const { unsupported("get(Index_t, std::vector<unsigned int> &)"); }
      virtual void get(Index_t, std::vector<unsigned long> &) const { unsupported("get(Index_t, std::vector<unsigned long> &)"); }
      virtual void get(Index_t, std::vector<std::string> &) const { unsupported("get(Index_t, std::vector<std::string> &)"); }

      virtual void set(Index_t, const bool &) { unsupported("get(Index_t, bool &)"); }
      virtual void set(Index_t, const double &) { unsupported("set(Index_t, const double &)"); }
      virtual void set(Index_t, const float &) { unsupported("set(Index_t, const float &)"); }
      virtual void set(Index_t, const char &) { unsupported("set(Index_t, const char &)"); }
      virtual void set(Index_t, const signed char &) { unsupported("set(Index_t, const signed char &)"); }
      virtual void set(Index_t, const signed short &) { unsupported("set(Index_t, const signed short &)"); }
      virtual void set(Index_t, const signed int &) { unsupported("set(Index_t, const signed int &)"); }
      virtual void set(Index_t, const signed long &) { unsupported("set(Index_t, const signed long &)"); }
      virtual void set(Index_t, const unsigned char &) { unsupported("set(Index_t, const unsigned char &)"); }
      virtual void set(Index_t, const unsigned short &) { unsupported("set(Index_t, const unsigned short &)"); }
      virtual void set(Index_t, const unsigned int &) { unsupported("set(Index_t, const unsigned int &)"); }
      virtual void set(Index_t, const unsigned long &) { unsupported("set(Index_t, const unsigned long &)"); }
      virtual void set(Index_t, const char *) { unsupported("set(Index_t, const char *)"); }
      virtual void set(Index_t, const std::string &) { unsupported("set(Index_t, const std::string &)"); }

      virtual void set(Index_t, const std::vector<bool> &) { unsupported("set(Index_t, const std::vector<bool> &)"); }
      virtual void set(Index_t, const std::vector<double> &) { unsupported("set(Index_t, const std::vector<double> &)"); }
      virtual void set(Index_t, const std::vector<float> &) { unsupported("set(Index_t, const std::vector<float> &)"); }
      virtual void set(Index_t, const std::vector<char> &) { unsupported("set(Index_t, const std::vector<char> &)"); }
      virtual void set(Index_t, const std::vector<signed char> &)
        { unsupported("set(Index_t, const std::vector<signed char> &)"); }
      virtual void set(Index_t, const std::vector<signed short> &)
        { unsupported("set(Index_t, const std::vector<signed short> &)"); }
      virtual void set(Index_t, const std::vector<signed int> &)
        { unsupported("set(Index_t, const std::vector<signed int> &)"); }
      virtual void set(Index_t, const std::vector<signed long> &)
        { unsupported("set(Index_t, const std::vector<signed long> &)"); }
      virtual void set(Index_t, const std::vector<unsigned char> &)
        { unsupported("set(Index_t, const std::vector<unsigned char> &)"); }
      virtual void set(Index_t, const std::vector<unsigned short> &)
        { unsupported("set(Index_t, const std::vector<unsigned short> &)"); }
      virtual void set(Index_t, const std::vector<unsigned int> &)
        { unsupported("set(Index_t, const std::vector<unsigned int> &)"); }
      virtual void set(Index_t, const std::vector<unsigned long> &)
        { unsupported("set(Index_t, const std::vector<unsigned long> &)"); }
      virtual void set(Index_t, const std::vector<std::string> &)
        { unsupported("set(Index_t, const std::vector<std::string> &)"); }

      virtual bool isNull(Index_t) const { unsupported("isNull() const"); return true; }
      virtual bool getNull(Index_t, bool &) const { unsupported("getNull(Index_t, bool &) const"); return true; }
      virtual bool getNull(Index_t, std::vector<bool> &) const
        { unsupported("getNull(Index_t, std::vector<bool> &) const"); return true; }

      /** \brief Copy a cell from another column to this column. An implementation which just throws an exception is
          provided in the base class to simplify derived classes which need only a small number of the methods to function.
      */
      virtual void copy(const IColumn *, Index_t, Index_t) { unsupported("copy(const IColumn *, Index_t, Index_t)"); }

      /** \brief Return a flag indicating whether this column holds scalar data.
      */
      virtual bool isScalar() const { return true; }

      /** \brief Return the name of the particular column implementation (subclass identifier).
      */
      virtual const std::string implementation() const = 0;

      /** \brief Get number of elements in the given cell. Default implementation assumes
          that only scalar valued cells are supported.
      */
      virtual Index_t getNumElements(Index_t = 0) const { return 1; }

      /** \brief Set number of elements in the given cell.
      */
      virtual void setNumElements(Index_t) { unsupported("setNumElements(Index_t)"); }

      /** \brief Get a string which identifies this column.
      */
      virtual const std::string & getId() const { return m_id; }

      /** \brief Get the units of this column, as a string.
      */
      virtual const std::string & getUnits() const { return m_units; }

      /** \brief Get a modifiable keyword associated with this column.
          \param base_name The base name of the keyword, which will be specialized for this column.
      */
      virtual Keyword & getColumnKeyword(const std::string &)
        { unsupported("getColumnKeyword"); return Keyword::emptyKeyword(); }

      /** \brief Get a constant keyword associated with this column.
          \param base_name The base name of the keyword, which will be specialized for this column.
      */
      virtual const Keyword & getColumnKeyword(const std::string &) const
        { unsupported("getColumnKeyword"); return Keyword::emptyKeyword(); }

      /// \brief Return a string identifying the full data type of the column.
      virtual std::string getFormat() const { unsupported("getFormat"); return ""; }

    protected:
      std::string m_units;

    private:
      /** \brief Private helper class to simplify throwing exceptions for unsupported features.
          \param method The name of the method throwing the exception.
      */
      void unsupported(const std::string & method) const {
        throw TipException(std::string("IColumn method ") + method + " is not supported for the " + implementation() +
          " implementation");
      }

      std::string m_id;
  };

}

#endif
