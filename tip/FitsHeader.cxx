/** \file FitsHeader.cxx

    \brief Implementation of utilities to help manage FITS specific table access.

    \author James Peachey, HEASARC
*/
#include <cctype>
#include <sstream>

#include "FitsHeader.h"
#include "tip/TipException.h"

namespace tip {

  FitsHeader::FitsHeader(const std::string & file_name, const std::string & ext_name,
    const std::string & filter, bool read_only): m_keyword_seq(), m_file_name(file_name), m_ext_name(ext_name),
    m_filter(filter), m_fp(0), m_is_primary(false), m_is_table(false), m_read_only(read_only) { open(); }

  // Close file automatically while destructing.
  FitsHeader::~FitsHeader() { close(); }

  // Subclasses call this to open the file and position it to the desired extension.
  void FitsHeader::open() {
    if (0 == m_fp) {
      fitsfile * fp = 0;
      int status = 0;

      // Construct the full name of the file from file name [extension] [filter] (if any):
      std::ostringstream s;
      s << m_file_name;
      if (!m_ext_name.empty()) s << "[" << m_ext_name << "]";
      if (!m_filter.empty()) s << "[" << m_filter << "]";
      std::string file_name = s.str();

      // Try to open the fits file read-write, unless read-only mode was explicitly set before open
      // was called.
      if (!m_read_only)
        fits_open_file(&fp, const_cast<char *>(file_name.c_str()), READWRITE, &status);

      // If opening read-write didn't work, or if read-only mode was explicitly set before open
      // was called...
      if (0 != status || m_read_only) {
        // Attempt to open the file read-only:
        status = 0;
        fits_open_file(&fp, const_cast<char *>(file_name.c_str()), READONLY, &status);
        m_read_only = true;
      }

      if (0 != status) {
        // TODO 9. 4/2/2004: Bug in cfitsio 2.48: Check for it and warn about it. The bug causes
        // the parser not to move to the correct extension.
        // DONE 9: 4/21/2004: A warning is issued for this specific version. Nothing more can be done, really.
        float cfitsio_version = 0.;
        fits_get_version(&cfitsio_version);
        // This is surreal. A FLOATING POINT VERSION NUMBER! Checking for == doesn't work -- I tried it.
        if (2.4795 < cfitsio_version && 2.4805 > cfitsio_version)
          throw TipException(status, std::string("WARNING: there is a known bug in Cfitsio 2.48's extended "
            "syntax parser!\nCould not open FITS extension ") + file_name);
        throw TipException(status, std::string("Could not open FITS extension \"") + file_name + '"');
      }

      // Success: save the pointer.
      m_fp = fp;

      // Read all keywords.
      loadAllKeywords();

      // See whether this is the primary extension.
      int hdu_num = 0;
      fits_get_hdu_num(m_fp, &hdu_num);
      m_is_primary = (1 == hdu_num);

      // Get the name of the extension. Do it this way, and not just by relying on m_ext_name in case
      // the user specified a wildcard.
      if (!m_is_primary) {
        try {
          getKeyword("EXTNAME", m_ext_name);
        } catch (const TipException &) {
        }
      }
      if (m_ext_name.empty()) {
        try {
          getKeyword("HDUNAME", m_ext_name);
        } catch (const TipException &) {
        }
      }

      // Check whether the file pointer is pointing at a table:
      int hdu_type = 0;
      fits_get_hdu_type(m_fp, &hdu_type, &status);
      if (0 != status) {
        close(status);
        throw TipException(status, formatWhat("Could not determine the type of the HDU"));
      }
      if (ASCII_TBL == hdu_type || BINARY_TBL == hdu_type)
        m_is_table = true;
    }
  }

  // Close file.
  void FitsHeader::close(int status) {
    if (0 != m_fp) {
      if (!m_read_only) fits_write_chksum(m_fp, &status);
      fits_close_file(m_fp, &status);
    }
    m_fp = 0;
  }

  Header::Iterator FitsHeader::find(const std::string & key_name) {
    Header::Iterator found_key = m_keyword_seq.begin();
    for (; found_key != m_keyword_seq.end() && key_name != found_key->getName(); ++found_key);
    return found_key;
  }

  Header::ConstIterator FitsHeader::find(const std::string & key_name) const {
    Header::ConstIterator found_key = m_keyword_seq.begin();
    for (; found_key != m_keyword_seq.end() && key_name != found_key->getName(); ++found_key);
    return found_key;
  }

  Header::Iterator FitsHeader::insert(Iterator itor, const KeyRecord & record) {
    int status = 0;
    fits_insert_record(m_fp, itor - m_keyword_seq.begin() + 1, const_cast<char *>(record.get().c_str()), &status);
    if (0 != status) {
      std::string msg = "Cannot insert record " + record.get();
      if (!itor->getName().empty()) msg += " before keyword " + itor->getName();
      throw TipException(status, formatWhat(msg));
    }
    return m_keyword_seq.insert(itor, record);
  }

  Header::Iterator FitsHeader::append(const KeyRecord & record) {
    return insert(m_keyword_seq.end(), record);
  }

  Header::Iterator FitsHeader::erase(Iterator itor) {
    int status = 0;
    fits_delete_record(m_fp, itor - begin() + 1, &status);
    return m_keyword_seq.erase(itor);
  }

  void FitsHeader::erase(const std::string & key_name) {
    int status = 0;
    // First, erase all matching keywords as far as cfitsio is concerned.
    do {
      fits_delete_key(m_fp, const_cast<char *>(key_name.c_str()), &status);
    } while (0 == status);
    if (KEY_NO_EXIST != status) throw TipException(status, formatWhat("Error deleting keyword \"" + key_name + "\""));

    // Next, erase all matching keywords in the container of keywords.
    // Iterate through container in reverse order, but using forward iterator.
    for (KeySeq_t::iterator itor = m_keyword_seq.end(); itor != m_keyword_seq.begin();) {
      // See if the previous item matches.
      if (key_name == (itor - 1)->getName()) {
        // Remove the matching item, and reset iterator to point to the element after the one
        // that was removed (end is OK because the iterator will be decremented above).
        itor = m_keyword_seq.erase(itor - 1);
      } else {
        --itor;
      }
    }
  }

  std::string FitsHeader::getKeyComment(const std::string & name) const {
    int status = 0;
    char value[FLEN_VALUE];
    char comment[FLEN_COMMENT];
    fits_read_keyword(m_fp, const_cast<char *>(name.c_str()), value, comment, &status);
    if (0 != status) throw TipException(status, formatWhat(std::string("Cannot read comment for keyword \"") + name + '"'));
    return comment;
  }

  void FitsHeader::setKeyComment(const std::string & name, const std::string & comment) {
    if (m_read_only)
      throw TipException(formatWhat(std::string("Cannot write comment for keyword \"") + name + "\"; object is not writable"));
    int status = 0;
    fits_modify_comment(m_fp, const_cast<char *>(name.c_str()), const_cast<char *>(comment.c_str()), &status);
    if (0 != status) throw TipException(status, formatWhat(std::string("Cannot write comment for keyword \"") + name + '"'));
  }

  std::string FitsHeader::getKeyUnit(const std::string & name) const {
    int status = 0;
    char unit[FLEN_CARD] = "";
    fits_read_key_unit(m_fp, const_cast<char *>(name.c_str()), unit, &status);
    if (0 != status) throw TipException(status, formatWhat(std::string("Cannot read unit for keyword \"") + name + '"'));
    return unit;
  }

  void FitsHeader::setKeyUnit(const std::string & name, const std::string & unit) {
    if (m_read_only)
      throw TipException(formatWhat(std::string("Cannot write unit for keyword \"") + name + "\"; object is not writable"));
    int status = 0;
    fits_write_key_unit(m_fp, const_cast<char *>(name.c_str()), const_cast<char *>(unit.c_str()), &status);
    if (0 != status) throw TipException(status, formatWhat(std::string("Cannot write unit for keyword \"") + name + '"'));
  }

  void FitsHeader::addComment(const std::string & comment) {
    if (m_read_only)
      throw TipException(formatWhat("Cannot add comment string; object is not writable"));
    int status = 0;
    fits_write_comment(m_fp, const_cast<char *>(comment.c_str()), &status);
    if (0 != status) throw TipException(status, formatWhat("Cannot add comment string"));
  }

  void FitsHeader::addHistory(const std::string & history) {
    if (m_read_only)
      throw TipException(formatWhat("Cannot add history string; object is not writable"));
    int status = 0;
    fits_write_history(m_fp, const_cast<char *>(history.c_str()), &status);
    if (0 != status) throw TipException(status, formatWhat("Cannot add history string"));
  }

  const std::string & FitsHeader::getName() const { return m_ext_name; }

  void FitsHeader::setName(const std::string & name) {
    if (m_read_only)
      throw TipException(formatWhat("Cannot set extension name; object is not writable"));
    if (m_is_primary) setKeyword("HDUNAME", name); else setKeyword("EXTNAME", name);
    m_ext_name = name;
  }

  std::string FitsHeader::formatWhat(const std::string & msg) const {
    std::ostringstream msg_str;
    msg_str << msg;
    if (!m_ext_name.empty()) msg_str << " in extension \"" << m_ext_name << '"';
    msg_str << " in file \"" << m_file_name << '"';
    return msg_str.str();
  }

  void FitsHeader::loadAllKeywords() {
    // Get number of keywords in header.
    int status = 0;
    int num_keywords = 0;
    fits_get_hdrspace(m_fp, &num_keywords, 0, &status);
    if (0 != status) throw TipException(status, formatWhat("Cannot get number of keywords in header"));

    // See if any of these keywords were not already loaded.
    if (num_keywords > 0 && (unsigned int)(num_keywords) > m_keyword_seq.size()) {
      int begin = m_keyword_seq.size();
      m_keyword_seq.resize(num_keywords);
      for (int ii = begin; ii < num_keywords; ++ii) {
        char card[FLEN_CARD] = "";
        // Read each keyword.
        fits_read_record(m_fp, ii + 1, card, &status);
        if (0 != status) {
          std::ostringstream os;
          os << "Cannot get keyword number " << ii;
          throw TipException(status, formatWhat(os.str()));
        }
        // Add key record to sequence.
        m_keyword_seq[ii] = KeyRecord(card);
      }
    }
  }

}
