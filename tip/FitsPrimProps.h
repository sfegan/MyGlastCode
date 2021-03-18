/** \file FitsPrimProps.h

    \brief Utilities to facilitate calling cfitsio functions from templated functions.

    \author James Peachey, HEASARC
*/
#ifndef tip_FitsPrimProps_h
#define tip_FitsPrimProps_h

namespace tip {

  /** \struct FitsPrimProps
      \brief Templated utility class holding type-specific information to simplify calling cfitsio functions
      from templated functions. This class is not part of the API and should not be of interest to clients.
  */
  template <typename T>
  struct FitsPrimProps {
      /** \brief Return the cfitsio code (TLONG, etc.) for the parameterized type.

          Note this is only defined for supported cfitsio primitives. An undefined symbol will result
          for types not already supported.
      */
      static int dataTypeCode();

      /** \brief Return the value used for undefined values in FITS files for the given primitive type.

          Note that this is tied to the primitive type being written, not to the type in the FITS file.
      */
      static T & undefined();
  };

}

#endif
