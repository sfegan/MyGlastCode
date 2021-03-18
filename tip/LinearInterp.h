/** \file LinearInterp.h

    \brief Utility to interpolate values in tables. Table must be ordered on the field used to interpolate.

    \author James Peachey, HEASARC
*/
#ifndef tip_LinearInterp_h
#define tip_LinearInterp_h

#include <string>
#include <vector>

#include "tip/Table.h"

namespace tip {

  /** \class LinearInterp

      \brief Utility to interpolate values in tables. Table must be ordered on the field used to interpolate.
  */
  class LinearInterp {
    public:
      LinearInterp(const Table::ConstIterator & begin, const Table::ConstIterator & end);

      void interpolate(const std::string & field, double value);

      double get(const std::string & field) const;

      void get(const std::string & field, std::vector<double> & value) const;

    private:
      Table::ConstIterator m_begin;
      Table::ConstIterator m_end;
      Table::ConstIterator m_match[2];
      double m_coeff[2];
  };

}

#endif
