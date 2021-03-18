//-*-mode:c++; mode:font-lock;-*-

/*! \file RandomNumbers.hpp

  RandomNumber selector - allow compile time use of the old supplier or
  the new generator

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 1.13 $
  \date       10/31/2007

  $Id: RandomNumbers.hpp 5429 2013-06-27 13:40:55Z sfegan $

*/

#ifndef RANDOMNUMBERS_HPP
#define RANDOMNUMBERS_HPP

#ifdef USE_RANDOMNUMBERS_TOS
#include "RandomNumbers_TOS.hpp"
typedef RandomNumbers_TOS RandomNumbers;
#else
#include "RandomNumbers_TNG.hpp"
typedef RandomNumbers_TNG<RNGCore::NR3Ran> RandomNumbers;
typedef SimpleRNGAdapter_TNG<RandomNumbers> SimpleRNGAdapter;
#endif

#endif
