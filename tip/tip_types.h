/** \file tip_types.h

    \brief Forward declarations and typedefs for all tip classes.

    \author James Peachey, HEASARC
*/
#ifndef tip_tip_types_h
#define tip_tip_types_h

namespace tip {

  /** \brief Type used to identify field index.
  */
  typedef signed int FieldIndex_t;

#ifdef TIP_USE_LONG_LONG_INDEX
  /** \brief Type used for differences in table row numbers.
  */
  typedef signed long long IndexDiff_t;

  /** \brief Type used for table row numbers.
  */
  typedef signed long long Index_t;

#else
  /** \brief Type used for differences in table row numbers.
  */
  typedef signed long IndexDiff_t;

  /** \brief Type used for table row numbers.
  */
  typedef signed long Index_t;
#endif

  /** \brief Type used to identify pixel ordinates.
  */
  typedef signed long PixOrd_t;
}

#endif
