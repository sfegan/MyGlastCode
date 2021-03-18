//-*-mode:c++; mode:font-lock;-*-

/*! \file SimpleRNG.hpp

  A simple RNG base class

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@llr.in2p3.fr         \n

  \version    $Revision: 1.13 $
  \date       09/24/2012

  $Id: SimpleRNG.hpp 5429 2013-06-27 13:40:55Z sfegan $

*/

#ifndef SIMPLERNG_HPP
#define SIMPLERNG_HPP

class SimpleRNG
{
public:
  virtual ~SimpleRNG();
  virtual double uniform() = 0;
};

#endif
