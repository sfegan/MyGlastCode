/** \file FitsColumn.h
    \brief FITS-specific interface to columnar data.
    \author James Peachey, HEASARC
*/
#ifndef tip_FitsColumn_h
#define tip_FitsColumn_h

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <memory>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "fitsio.h"

#include "FitsTable.h"
#include "FitsPrimProps.h"
#include "tip/IColumn.h"
#include "tip/TipException.h"
#include "tip/tip_types.h"

namespace tip {

  template <typename T>
  class FitsColumn : public IColumn {
    public:
      FitsColumn(FitsTable * ext, const std::string & id, FieldIndex_t field_index);

      virtual ~FitsColumn() throw() {}

      virtual void get(Index_t record_index, double & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, float & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, char & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed char & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed short & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed int & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, signed long & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned char & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned short & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned int & dest) const { getScalar(record_index, dest); }
      virtual void get(Index_t record_index, unsigned long & dest) const { getScalar(record_index, dest); }

      virtual void get(Index_t record_index, std::vector<double> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<float> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<char> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed char> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed short> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed int> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<signed long> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned char> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned short> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned int> & dest) const { getVector(record_index, dest); }
      virtual void get(Index_t record_index, std::vector<unsigned long> & dest) const { getVector(record_index, dest); }

      virtual void set(Index_t record_index, const double & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const float & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const char & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed char & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed short & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed int & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const signed long & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned char & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned short & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned int & src) { setScalar(record_index, src); }
      virtual void set(Index_t record_index, const unsigned long & src) { setScalar(record_index, src); }

      virtual void set(Index_t record_index, const std::vector<double> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<float> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<char> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed char> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed short> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed int> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<signed long> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned char> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned short> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned int> & src) { setVector(record_index, src); }
      virtual void set(Index_t record_index, const std::vector<unsigned long> & src) { setVector(record_index, src); }

      virtual void get(Index_t record_index, std::string & dest) const {
        char * buf = 0;
        try {
          // For strings, make a buffer to hold the value.
          buf = new char[m_display_width + 1];

          // Set the buffer contents to 0.
          std::memset(buf, '\0', m_display_width + 1);

          // Now call fitsio to fill the buffer with the string.
          getScalar(record_index, buf);

          // Copy this to the input.
          dest = buf;
        } catch (...) {
          delete [] buf;
          throw;
        }
        delete [] buf;
      }

      virtual void get(Index_t record_index, std::vector<std::string> & dest) const {
        // Clear content of destination in case there is a problem.
        dest.clear();
        char * buf(new char[m_display_width + 1]);
        char * tmp_dest[] = { buf };

        // Redimension the destination so that it will hold the contents of this vector.
        dest.resize(m_repeat);
        
        for (std::vector<std::string>::size_type index = 0; index != dest.size(); ++index) {
          // Read each element in the column, letting cfitsio do the conversions to strings.
          int status = 0;
          int any_null = 0;
          fits_read_col(m_ext->getFp(), FitsPrimProps<char *>::dataTypeCode(), m_field_index, record_index + 1,
            index + 1, 1, 0, tmp_dest, &any_null, &status);
          if (0 != status) {
            std::ostringstream os;
            os << "FitsColumn::get(Index_t, std::vector<std::string> &) could not read record " << record_index;
            delete [] buf;
            throw TipException(status, os.str());
          }
          if (0 != any_null) {
            dest[index] = FitsPrimProps<char *>::undefined();
          } else {
            // Skip leading and trailing spaces.
            char * begin = *tmp_dest;
            while ('\0' != *begin && 0 != std::isspace(*begin)) ++begin;
            char * end = begin;
            while ('\0' != *end && 0 == std::isspace(*end)) ++end;

            // Copy cfitsio's string to the output string.
            dest[index].assign(begin, end);
          }
        }
        delete [] buf;
      }

      virtual void set(Index_t record_index, const char * src) { set(record_index, std::string(src)); }

      virtual void set(Index_t record_index, const std::string & src) {
        // Writing strings must be handled carefully, because cfitsio only converts strings to numbers when reading.
        switch (m_type_code) {
          case TSTRING:
            // String data -> string column.
            setScalar(record_index, src.c_str());
            break;
          case TLOGICAL: {
              // String data -> boolean column.
              bool logical_val = false;
              std::string src_copy(src);
              for (std::string::iterator itor = src_copy.begin(); itor != src_copy.end(); ++itor) *itor = tolower(*itor);
              if (src_copy == "t" || src_copy == "true") logical_val = true;
              else if (src_copy == "f" || src_copy == "false") logical_val = false;
              else throw TipException(
                "FitsColumn::set(Index_t, const std::string &) could not convert string \"" + src + "\" to bool");
              set(record_index, logical_val);
              break;
            }
          default: {
              // String data -> numeric column. Use double for everything.
              char * remainder = 0;
              double dval = strtod(src.c_str(), &remainder);
              if (0 != remainder && '\0' != *remainder) throw TipException(
                "FitsColumn::set(Index_t, const std::string &) could not convert string \"" + src + "\" to double");
              set(record_index, dval);
              break;
            }
        }
      }

      virtual void set(Index_t record_index, const std::vector<std::string> & src) {
        if (m_scalar) throw TipException("FitsColumn::set(Index_t, const vector<string> &) was called but field is not a vector");
        // Writing strings must be handled carefully, because cfitsio only converts strings to numbers when reading.
        switch (m_type_code) {
          case TSTRING: {
              // String data -> string column. Copy pointers to individual strings in src vector into an array of C primitives
              // to pass to cfitsio.
              std::vector<const char *> buf(src.size());
              for (std::vector<std::string>::size_type idx = 0; idx != src.size(); ++idx) {
                buf[idx] = &src[idx][0];
              }
              setVector(record_index, buf);
              break;
            }
          case TLOGICAL: {
              // String data -> boolean column.
              std::vector<char> buf(src.size());
              for (std::vector<std::string>::size_type idx = 0; idx != src.size(); ++idx) {
                // Strip leading/trailing blanks.
                std::string::size_type begin = src[idx].find_first_not_of(" \n\t\v");
                std::string::size_type end = src[idx].find_last_of(" \n\t\v");
                if (begin == std::string::npos) begin = 0;
                if (end == std::string::npos) end = src[idx].size();
                std::string src_copy(src[idx].begin() + begin, src[idx].begin() + end);
                for (std::string::iterator itor = src_copy.begin(); itor != src_copy.end(); ++itor) *itor = tolower(*itor);
                if (src_copy == "t" || src_copy == "true") buf[idx] = true;
                else if (src_copy == "f" || src_copy == "false") buf[idx] = false;
                else if (src_copy.empty()) buf[idx] = FitsPrimProps<char>::undefined();
                else throw TipException(
                  "FitsColumn::set(Index_t, const vector<string> &) could not convert string \"" + src[idx] + "\" to bool");
              }
              setVector(record_index, buf);
              break;
            }
          default: {
              // String data -> numeric column. Use double for everything.
              std::vector<double> buf(src.size());
              char * remainder = 0;
              for (std::vector<std::string>::size_type idx = 0; idx != src.size(); ++idx) {
                buf[idx] = strtod(src[idx].c_str(), &remainder);
                if (0 != remainder && '\0' != *remainder) throw TipException(
                  "FitsColumn::set(Index_t, const vector<string> &) could not convert string \"" + src[idx] + "\" to double");
              }
              setVector(record_index, buf);
              break;
            }
        }
      }

      // Specializations for bool. This is done instead of specializing the templates to avoid complexity and
      // platform-specific code.
      virtual void get(Index_t record_index, bool & dest) const {
        if (!m_scalar) throw TipException("FitsColumn::get(Index_t, bool &) was called but field is not a scalar");
        int status = 0;
        char tmp_dest = 0;
        int any_nul = 0;
        fits_read_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          0, &tmp_dest, &any_nul, &status);
        if (0 != status) throw TipException(status, "FitsColumn::get(Index_t, bool &) failed to read scalar cell value");
        if (0 != any_nul) dest = FitsPrimProps<bool>::undefined();
        dest = (0 != tmp_dest);
      }

      virtual void get(Index_t record_index, std::vector<bool> & dest) const {
        if (m_scalar) throw TipException("FitsColumn::get(Index_t, std::vector<bool> &) was called but field is not a vector");
        int status = 0;
        Index_t num_els = getNumElements(record_index);
        char * tmp_dest = new char[num_els];
        int * any_nul = new int[num_els];
        fits_read_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          0, tmp_dest, any_nul, &status);
        if (0 != status) {
          delete [] any_nul;
          delete [] tmp_dest;
          throw TipException(status, "FitsColumn::get(Index_t, std::vector<bool> &) failed to read vector cell value");
        }
        for (Index_t ii = 0; ii != num_els; ++ii)
          if (0 == any_nul[ii]) dest[ii] = (0 != tmp_dest[ii]); else dest[ii] = FitsPrimProps<bool>::undefined();
        delete [] any_nul;
        delete [] tmp_dest;
      }

      virtual void set(Index_t record_index, const bool & src) {
        if (!m_scalar) throw TipException("FitsColumn::set(Index_t, const bool &) called but field is not a scalar");
        int status = 0;
        char tmp_src = src;
        fits_write_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          const_cast<void *>(static_cast<const void *>(&tmp_src)), &status);
        if (0 != status) throw TipException(status, "FitsColumn::set(Index_t, const bool &) failed to write scalar cell value");
      }

      virtual void set(Index_t record_index, const std::vector<bool> & src) {
        if (m_scalar)
          throw TipException("FitsColumn::set(Index_t, const std::vector<bool> &) called but field is not a vector");
        int status = 0;
        Index_t num_els = src.size();

        if (!m_var_length && num_els > m_repeat) {
          std::ostringstream os;
          os << "FitsColumn::set(Index_t, const std::vector<bool> &) attempted to write " << num_els <<
            " elements into a cell of size " << m_repeat;
          throw TipException(os.str());
        }

        char * tmp_src = new char[num_els];
        for (Index_t ii = 0; ii < num_els; ++ii) tmp_src[ii] = src[ii];
        fits_write_col(m_ext->getFp(), FitsPrimProps<bool>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          tmp_src, &status);
        delete [] tmp_src;
        if (0 != status)
          throw TipException(status, "FitsColumn::set(Index_t, const std::vector<bool> &) failed to write vector cell value");
      }

      virtual bool isNull(Index_t record_index) const {
        if (!m_scalar) throw TipException("FitsColumn::isNull(Index_t) called but field is not a scalar");
        int status = 0;
        int any_null = 0;
        // For strings, make a buffer to hold the value.
        char * dest = new char[m_display_width + 1];
        fits_read_col(m_ext->getFp(), FitsPrimProps<char *>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          &FitsPrimProps<char *>::undefined(), &dest, &any_null, &status);
        delete [] dest; dest = 0;
        if (0 != status) throw TipException(status, "FitsColumn::isNull failed to read scalar cell value");
        return 0 != any_null;
      }

      virtual bool getNull(Index_t record_index, bool & null_value) const {
        return null_value = isNull(record_index);
      }

      virtual bool getNull(Index_t record_index, std::vector<bool> & null_value) const {
        if (m_scalar) throw TipException("FitsColumn::getNull(Index_t, std::vector<bool> &) called but field is not a vector");
        int status = 0;
        int any_null = 0;
        Index_t num_els = getNumElements(record_index);
        null_value = std::vector<bool>(num_els, false);
        // Array to hold null values in a form native to cfitsio.
        char * null_value_cp = new char[num_els];
        // For strings, make a buffer to hold the values.
        char * buf = new char[num_els * (m_display_width + 1)];
        // Make an destination array of pointers to the individual sub-buffers.
        char ** dest = new char *[num_els];
        // Set pointers in destination array.
        for (Index_t ii = 0; ii != num_els; ++ii) dest[ii] = buf + ii * (m_display_width + 1);

        // Read the column value, looking for nulls.
        fits_read_colnull(m_ext->getFp(), FitsPrimProps<char *>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          dest, null_value_cp, &any_null, &status);

        // Copy any nulls found to the output array.
        //if (any_null) for (Index_t ii = 0; ii != num_els; ++ii) null_value[ii] = 0 != null_value_cp[ii];
        for (Index_t ii = 0; ii != num_els; ++ii) null_value[ii] = 0 != null_value_cp[ii];

        // Clean up.
        delete [] dest; dest = 0;
        delete [] buf; buf = 0;
        delete [] null_value_cp; null_value_cp = 0;

        if (0 != status) throw TipException(status, "FitsColumn::getNull failed to read vector cell value");
        return 0 != any_null;
      }

      /** \brief Copy a cell from another column to this column.
          \param src Pointer to the source column.
          \param src_index Index of the cell in the source column.
          \param dest_index Index of the cell in the destination column (this column).
      */
      virtual void copy(const IColumn * src, Index_t src_index, Index_t dest_index) {
        if (m_scalar) {
          T val;
          src->get(src_index, val);
          set(dest_index, val);
        } else {
          std::vector<T> val;
          src->get(src_index, val);
          set(dest_index, val);
        }
      }

      /** \brief Return a flag indicating whether this column holds scalar data.
      */
      virtual bool isScalar() const { return m_scalar; }

      /** \brief Return the name of the particular column implementation.
      */
      virtual const std::string implementation() const  { return "FITS"; }

      /** \brief Get number of elements in the given cell.
          \param record_index The record number identifying the cell.
      */
      virtual Index_t getNumElements(Index_t record_index = 0) const {
        if (!m_var_length) return m_repeat;
        int status = 0;
        Index_t num_els = 0;
        // Get number of elements in this particular field.
#ifdef TIP_USE_LONG_LONG_INDEX
        fits_read_descriptll(m_ext->getFp(), m_field_index, record_index + 1, &num_els, 0, &status);
#else
        fits_read_descript(m_ext->getFp(), m_field_index, record_index + 1, &num_els, 0, &status);
#endif
        if (0 != status) throw TipException(status, "FitsColumn::getNumElements failed to get size of variable length cell");

        return num_els;
      }

      /** \brief Set number of elements in the given cell.
      */
      virtual void setNumElements(Index_t num_elements) {
        if (m_var_length) throw TipException("FitsColumn::setNumElements cannot change the width of variable length column");
        int status = 0;
        fits_modify_vector_len(m_ext->getFp(), m_field_index, num_elements, &status);
        if (0 != status) throw TipException(status, "FitsColumn::setNumElements failed to modify field");

        // Update column information.
        m_repeat = num_elements;
        if (1 == m_repeat && !m_var_length) m_scalar = true; else m_scalar = false;
      }

      /** \brief Get a modifiable keyword associated with this column.
          \param base_name The base name of the keyword, which will be specialized for this column.
      */
      virtual Keyword & getColumnKeyword(const std::string & base_name);

      /** \brief Get a constant keyword associated with this column.
          \param base_name The base name of the keyword, which will be specialized for this column.
      */
      virtual const Keyword & getColumnKeyword(const std::string & base_name) const;

      /// \brief Return a string identifying the full data type of the column.
      virtual std::string getFormat() const { return m_type_string; }

    private:
      template <typename U>
      void getScalar(Index_t record_index, U & dest) const {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (!m_scalar) throw TipException("FitsColumn::getScalar was called but field is not a scalar");
        int status = 0;
        int any_null = 0;
        fits_read_col(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          &FitsPrimProps<U>::undefined(), &dest, &any_null, &status);
        if (0 != status) throw TipException(status, "FitsColumn::getScalar failed to read scalar cell value");
      }

      template <typename U>
      void getVector(Index_t record_index, std::vector<U> & dest) const {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (m_scalar) throw TipException("FitsColumn::getVector was called but field is not a vector");
        int status = 0;
        Index_t num_els = getNumElements(record_index);
        dest.resize(num_els);
        U * dest_begin = &dest.front();
        int any_null = 0;
        fits_read_col(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          &FitsPrimProps<U>::undefined(), dest_begin, &any_null, &status);
        if (0 != status) throw TipException(status, "FitsColumn::getVector failed to read vector cell value");
      }

      template <typename U>
      void setScalar(Index_t record_index, const U & dest) {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (!m_scalar) throw TipException("FitsColumn::setScalar called but field is not a scalar");
        int status = 0;
        if (m_ext->readOnly()) throw TipException("FitsColumn::setScalar called for a read-only file");
        fits_write_colnull(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, m_repeat,
          const_cast<void *>(static_cast<const void *>(&dest)), &FitsPrimProps<U>::undefined(), &status);
        if (0 != status) throw TipException(status, "FitsColumn::setScalar failed to write scalar cell value");
      }

      template <typename U>
      void setVector(Index_t record_index, const std::vector<U> & src) {
        // Prevent accidental calling for bool or string. The optimizer will swallow this.
        assert(typeid(U) != typeid(bool) && typeid(U) != typeid(std::string));
        if (m_scalar) throw TipException("FitsColumn::setVector called but field is not a vector");
        if (m_ext->readOnly()) throw TipException("FitsColumn::setVector called for a read-only file");
        int status = 0;
        Index_t num_els = src.size();

        if (!m_var_length && num_els > m_repeat) {
          std::ostringstream os;
          os << "FitsColumn::setVector attempted to write " << num_els << " elements into a cell of size " << m_repeat;
          throw TipException(os.str());
        }

        const U * src_begin = &src.front();
        fits_write_colnull(m_ext->getFp(), FitsPrimProps<U>::dataTypeCode(), m_field_index, record_index + 1, 1, num_els,
          const_cast<void *>(static_cast<const void *>(src_begin)), &FitsPrimProps<U>::undefined(), &status);
        if (0 != status) throw TipException(status, "FitsColumn::setVector failed to write vector cell value");
      }

      std::string m_type_string;
      FitsTable * m_ext;
      FieldIndex_t m_field_index;
      Index_t m_repeat;
      Index_t m_width;
      int m_type_code;
      int m_display_width;
      bool m_var_length;
      bool m_scalar;
  };

  template <typename T>
  inline FitsColumn<T>::FitsColumn(FitsTable * ext, const std::string & id, FieldIndex_t field_index): IColumn(id),
    m_type_string(), m_ext(ext), m_field_index(field_index), m_repeat(0), m_width(0), m_type_code(0), m_display_width(0),
    m_var_length(false), m_scalar(false) {

    // Determine characteristics of this column.
    int status = 0;
#ifdef TIP_USE_LONG_LONG_INDEX
    fits_get_coltypell(m_ext->getFp(), m_field_index, &m_type_code, &m_repeat, &m_width, &status);
#else
    fits_get_coltype(m_ext->getFp(), m_field_index, &m_type_code, &m_repeat, &m_width, &status);
#endif
    if (0 != status) throw TipException(status, "FitsColumn::FitsColumn failed to get format of field " + id);

    // Read the TFORM keyword, which gives the full layout of the column.
    getColumnKeyword("TFORM").get(m_type_string);

    // Detect unsigned integral types, and modify m_type_string as needed to reflect this.
    if (TINT == m_type_code || TLONG == m_type_code || TSHORT == m_type_code) {
      // Check the tscal.
      double tscal = 0.;
      try {
        getColumnKeyword("TSCAL").get(tscal);
      } catch (const TipException &) {
        tscal = 1.;
      }

      // If tscal is 1., need to check TZERO to see if it indicates the correct offset for unsigned integers.
      if (1. == tscal) {
        unsigned long tzero = 0;
        try {
          getColumnKeyword("TZERO").get(tzero);
          // Check whether TZERO is the correct offset for representing unsigned integers.
          switch (m_type_code) {
            case TINT:
            case TLONG:
              if (1u<<31u == tzero) {
                std::string::size_type index = m_type_string.find("J");
                if (std::string::npos != index) m_type_string[index] = 'V';
              }
              break;
            case TSHORT:
              if (1u<<15u == tzero) {
                std::string::size_type index = m_type_string.find("I");
                if (std::string::npos != index) m_type_string[index] = 'U';
              }
              break;
            default:
              break;
          }
        } catch (const TipException &) {
          // Failure to read TZERO as an unsigned long implies only that this column/field does not contain
          // an unsigned integer, so just continue.
        }
      }
    }

    // Handle special case of strings, for which the info returned by fits_get_coltype means something different.
    if (TSTRING == m_type_code) {
      m_repeat /= m_width;
      if (0 == m_repeat) m_repeat = 1;
    }

    // Handle variable length columns.
    if (m_type_code < 0) {
      m_type_code *= -1;
      m_var_length = true;
    }

    // Determine whether this is a scalar column.
    if (1 == m_repeat && !m_var_length) m_scalar = true;

    // Get units.
    std::ostringstream os;
    os << "TUNIT" << m_field_index;
    char units[FLEN_CARD] = "";
    fits_read_key(m_ext->getFp(), TSTRING, const_cast<char *>(os.str().c_str()), units, 0, &status);
    if (0 == status) m_units = units;
    else if (KEY_NO_EXIST != status) throw TipException(status, "FitsColumn::FitsColumn failed to get units of field");

    // Find out how big strings need to be to accomodate this data field.
    status = 0;
    fits_get_col_display_width(m_ext->getFp(), m_field_index, &m_display_width, &status);
    if (0 != status) throw TipException(status, "FitsColumn::FitsColumn failed to get string width of field " + id);
  }

  template <typename T>
  inline Keyword & FitsColumn<T>::getColumnKeyword(const std::string & base_name) {
    std::ostringstream os;
    os << base_name << m_field_index;
    return m_ext->getHeader()[os.str()];
  }

  template <typename T>
  inline const Keyword & FitsColumn<T>::getColumnKeyword(const std::string & base_name) const {
    std::ostringstream os;
    os << base_name << m_field_index;
    return m_ext->getHeader()[os.str()];
  }

}

#endif
