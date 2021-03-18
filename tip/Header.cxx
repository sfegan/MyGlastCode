/** \file Header.cxx

    \brief High level encapsulation of a FITS-like header.

    \author James Peachey, HEASARC
*/
#include <exception>
#include "tip/Header.h"

namespace tip {

  Keyword & Keyword::emptyKeyword() {
    static Keyword s_empty_keyword;
    return s_empty_keyword;
  }

  void Header::get(const char ** keys, KeyValCont_t & key_vals) const {
    std::string value;
    // Iterate over the given set of keys.
    for (const char ** key = keys; 0 != *key; ++key) {
      try {
        // Look up its value.
        operator [](*key).get(value);

        // Put it into the output group.
        key_vals.push_back(KeyValPair_t(*key, value));

      } catch (const std::exception &) {
        // Ignore problem getting a keyword.
      }
    }
  }

  void Header::update(const KeyValCont_t & key_vals) {
    std::string value;
    for (KeyValCont_t::const_iterator itor = key_vals.begin(); itor != key_vals.end(); ++itor) {
      try {
        // Find the keyword.
        Keyword & keyword = operator [](itor->first);

        // Read the keyword. This will throw if the keyword does not already exist.
        keyword.get(value);

        // Update the keyword.
        keyword.set(itor->second);
      } catch (const std::exception &) {
        // Ignore problem getting a keyword, but do not set its value.
      }
    }
  }

  std::string Header::formatTime(const time_t & time) const {
    // Standard date format defined by FITS standard.
    char string_time[] = "YYYY-MM-DDThh:mm:ss";

    // Format using ctime functions.
    struct tm * loc_time = localtime(&time);
    strftime(string_time, sizeof(string_time), "%Y-%m-%dT%H:%M:%S", loc_time);

    // Return formatted time string.
    return string_time;
  }

}
