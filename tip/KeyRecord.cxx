/** \file KeyRecord.cxx
    \brief Implementation of KeyRecord class.
    \authors Lawrence Brown, HEASARC/GSSC
             James Peachey, HEASARC/GSSC
*/
#include <cctype>
#include <errno.h>
#include <cstdlib>
#include <sstream>
#include <cstring>

#include "fitsio.h"
#include "tip/KeyRecord.h"
#include "tip/TipException.h"

// Cfitsio doesn't define this, so we will:
#define FLEN_KEYNAME 9

namespace tip {

  KeyRecord::KeyRecord(const std::string & record): m_record(record) {}

  // Check whether a value is defined.
  bool KeyRecord::empty() const { return getValue().empty(); }

  const std::string & KeyRecord::get() const { return m_record; }

  void KeyRecord::set(const std::string & record) { m_record = record; }

  std::string KeyRecord::getName() const {
    char key_name[FLEN_CARD];
    int status = 0;
    int len = 0;
    fits_get_keyname(const_cast<char *>(m_record.c_str()), key_name, &len, &status);
    if (0 != status) throw TipException(status, "KeyRecord::getName could not parse record");
    return std::string(key_name, key_name + len);
  }

  std::string KeyRecord::getValue() const {
    char value[FLEN_VALUE] = "";
    int status = 0;
    fits_parse_value(const_cast<char *>(m_record.c_str()), value, 0, &status);
    if (0 != status) throw TipException(status, "KeyRecord::getValue could not parse record");

    // See if this is treated as a string, in which case it will have trailing blanks and a quote.
    char * ptr = value + strlen(value);

    // Skip trailing space outside quote.
    while(ptr > value && 0 != std::isspace(*(ptr - 1))) --ptr;

    // Handle the trailing quote if any.
    if (ptr > value && '\'' == *(ptr - 1)) {
      // Skip trailing quote.
      --ptr;

      // Skip trailing space inside quote.
      while(ptr > value && 0 != std::isspace(*(ptr - 1))) --ptr;
    }

    // Make sure string is properly terminated at this point.
    *ptr = '\0';

    ptr = value;

    // Handle leading quote, if any.
    if ('\'' == *ptr) ++ptr;

    return ptr;
  }

  void KeyRecord::setValue(const std::string & value) {
    char orig_value[FLEN_VALUE] = ""; // Key name length in fits <= 8.
    char name[FLEN_KEYNAME] = ""; // Key name length in fits <= 8.
    char comment[FLEN_COMMENT] = "";
    int status = 0;
    int key_length = 0;

    // Get name and comment.
    fits_get_keyname(const_cast<char *>(m_record.c_str()), name, &key_length, &status);
    if (0 != status) throw TipException(status, "KeyRecord::setValue could not get name of key record");
    else if (!isalnum(*name)) throw TipException("KeyRecord::setValue: name of key record starts with non-alpha-numeric character");
    fits_parse_value(const_cast<char *>(m_record.c_str()), orig_value, comment, &status);
    if (0 != status) throw TipException(status, "KeyRecord::setValue could not parse record");

    // Classifying state of this record.
    // Fits standard states a ' must be in column 11 to indicate a string.
    bool is_orig_string = false;
    if (m_record.size() > 10 && m_record.substr(10, 1) == "'") is_orig_string = true;
    bool is_orig_blank = ('\0' == *orig_value);
    bool is_orig_num = !is_orig_string && !is_orig_blank;

    // Peculiarity of cfitsio. If keyword is blank, the comment includes leading spaces followed by "/ ",
    // so remove them in this case only.
    char * com_ptr = comment;
    if (is_orig_blank) {
      // Remove leading spaces.
      while(isspace(*com_ptr)) ++com_ptr;
      // Remove / and the next character (assumed to be a space).
      if ('/' == *com_ptr) { ++com_ptr; ++com_ptr; }
    }

    // Classify the type of the new value.
    // First see if value looks like a logical (hence not a string).
    bool is_new_string = true;
    if (value == "T" || value == "F") is_new_string = false;
    else {
      // Does not look like a logical, so see if it looks like a number.
      char * remainder = 0;
      errno = 0;
      strtod(value.c_str(), &remainder);
      if (0 == errno && (0 == remainder || '\0' == *remainder)) is_new_string = false;
    }

    // Truth table to determine how to treat the new string.
    bool treat_as_string = true;
    if (is_orig_blank && !is_new_string) treat_as_string = false;
    if (is_orig_num && !is_new_string) treat_as_string = false;
    if (is_orig_num && is_new_string) throw TipException("KeyRecord::setValue attempted to assign a string to a numeric key record");

    // Start constructing the new record.
    std::string new_rec;

    std::ostringstream os;
    os.width(8);
    os << std::left << name << "= ";

    if(treat_as_string) {
      // Add leading quote.
      os << "'";

      // Add value, again left justified, padded as needed for fixed format FITS standard keywords.
      os.width(8);
      os << value << "'";

      // Save result so far.
      std::string result = os.str();
      os.str("");

      // Combine comment with its delimiter.
      std::string com = std::string(" / ") + com_ptr;

      // Want comment to be padded by up to 10 characters on the left, respecting overall length of record.
      os.width(50);
      os << std::left << com;
      com = os.str();
      os.str("");

      os.width(80 - result.size());
      os << std::right << com;
      new_rec = result + os.str();

    } else {
      // Not a string: very straightforward fixed length formatting.
      os.width(20);
      os << std::right << value << std::left << " / " << com_ptr;
      new_rec = os.str();
    }

    // Restrict overall length.
    new_rec = new_rec.substr(0, FLEN_CARD - 1);

    // Remove trailing spaces.
    int pos = new_rec.find_last_not_of(" ");
    new_rec.erase(pos + 1);

    // Finally replace the record.
    m_record = new_rec;
  }

  std::string KeyRecord::getComment() const {
    char value[FLEN_VALUE];
    char comment[FLEN_COMMENT];
    int status = 0;
    fits_parse_value(const_cast<char *>(m_record.c_str()), value, comment, &status);
    if (0 != status) throw TipException(status, "KeyRecord::getComment could not parse record");

    // In some circumstances Cfitsio returns leading spaces and/or leading / in the comment. Go figure.
    char * cp = comment;
    while(0 != isspace(*cp)) ++cp;
    if ('/' == *cp) ++cp;
    while(0 != isspace(*cp)) ++cp;
    return cp;
  }

}
