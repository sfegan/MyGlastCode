/** \file FitsPrimProps.cxx

    \brief Templated utility class holding type-specific information to simplify calling cfitsio functions
    from templated functions. This class is not part of the API and should not be of interest to clients.

    \author James Peachey, HEASARC
*/
#include "fitsio.h"

#include "FitsPrimProps.h"

#include <limits>
#include <string>

namespace {
  char * s_cp_undef = "INDEF";
  const char * s_ccp_undef = "INDEF";
  bool s_bool_undef = std::numeric_limits<bool>::min();
  double s_double_undef = std::numeric_limits<double>::min();
  float s_float_undef = std::numeric_limits<float>::min();
  char s_char_undef = std::numeric_limits<char>::min();
  signed char s_signed_char_undef = std::numeric_limits<signed char>::min();
  signed short s_signed_short_undef = std::numeric_limits<signed short>::min();
  signed int s_signed_int_undef = std::numeric_limits<signed int>::min();
  signed long s_signed_long_undef = std::numeric_limits<signed long>::min();
  unsigned char s_unsigned_char_undef = std::numeric_limits<unsigned char>::min();
  unsigned short s_unsigned_short_undef = std::numeric_limits<unsigned short>::min();
  unsigned int s_unsigned_int_undef = std::numeric_limits<unsigned int>::min();
  unsigned long s_unsigned_long_undef = std::numeric_limits<unsigned long>::min();
  long long s_long_long_undef = std::numeric_limits<long long>::min();
}

namespace tip {

  // Currently these codes are not defined for complex and double complex cfitsio types.
  template <> int FitsPrimProps<char *>::dataTypeCode() { return TSTRING; }
  template <> int FitsPrimProps<const char *>::dataTypeCode() { return TSTRING; }
  template <> int FitsPrimProps<std::string>::dataTypeCode() { return TSTRING; }
  template <> int FitsPrimProps<bool>::dataTypeCode() { return TLOGICAL; }
  template <> int FitsPrimProps<double>::dataTypeCode() { return TDOUBLE; }
  template <> int FitsPrimProps<float>::dataTypeCode() { return TFLOAT; }
  template <> int FitsPrimProps<char>::dataTypeCode() { return TBYTE; }
  template <> int FitsPrimProps<signed char>::dataTypeCode() { return TBYTE; }
  template <> int FitsPrimProps<signed short>::dataTypeCode() { return TSHORT; }
  template <> int FitsPrimProps<signed int>::dataTypeCode() { return TINT; }
  template <> int FitsPrimProps<signed long>::dataTypeCode() { return TLONG; }
  template <> int FitsPrimProps<unsigned char>::dataTypeCode() { return TBYTE; }
  template <> int FitsPrimProps<unsigned short>::dataTypeCode() { return TUSHORT; }
  template <> int FitsPrimProps<unsigned int>::dataTypeCode() { return TUINT; }
  template <> int FitsPrimProps<unsigned long>::dataTypeCode() { return TULONG; }
  template <> int FitsPrimProps<long long>::dataTypeCode() { return TLONGLONG; }

  template <> char * & FitsPrimProps<char *>::undefined() { return s_cp_undef; }
  template <> const char * & FitsPrimProps<const char *>::undefined() { return s_ccp_undef; }
  template <> bool & FitsPrimProps<bool>::undefined() { return s_bool_undef; }
  template <> double & FitsPrimProps<double>::undefined() { return s_double_undef; }
  template <> float & FitsPrimProps<float>::undefined() { return s_float_undef; }
  template <> char & FitsPrimProps<char>::undefined() { return s_char_undef; }
  template <> signed char & FitsPrimProps<signed char>::undefined() { return s_signed_char_undef; }
  template <> signed short & FitsPrimProps<signed short>::undefined() { return s_signed_short_undef; }
  template <> signed int & FitsPrimProps<signed int>::undefined() { return s_signed_int_undef; }
  template <> signed long & FitsPrimProps<signed long>::undefined() { return s_signed_long_undef; } 
  template <> unsigned char & FitsPrimProps<unsigned char>::undefined() { return s_unsigned_char_undef; }
  template <> unsigned short & FitsPrimProps<unsigned short>::undefined() { return s_unsigned_short_undef; }
  template <> unsigned int & FitsPrimProps<unsigned int>::undefined() { return s_unsigned_int_undef; }
  template <> unsigned long & FitsPrimProps<unsigned long>::undefined() { return s_unsigned_long_undef; }
  template <> long long & FitsPrimProps<long long>::undefined() { return s_long_long_undef; }

}
