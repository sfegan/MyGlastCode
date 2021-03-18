// RandomNumbers.i - SWIG wrapper for the RandomNumbers classes
// Stephen Fegan - sfegan@llr.in2p3.fr - June 2013
// $Id: RandomNumbers.i 5429 2013-06-27 13:40:55Z sfegan $

%module RandomNumbers
%{
#include "SimpleRNG.hpp"
#include "RandomNumbers.hpp"
#include "RandomNumbers_TNG.hpp"
%}

%include "SimpleRNG.hpp"
%include "RandomNumbers_TNG.hpp"
%template(RandomNumbers) RandomNumbers_TNG<RNGCore::NR3Ran>;
%template(SimpleRNGAdapter) SimpleRNGAdapter_TNG<RandomNumbers>;
%include "RandomNumbers.hpp"
