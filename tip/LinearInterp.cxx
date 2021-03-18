/** \file LinearInterp.cxx

    \brief Utility to interpolate values in tables. Table must be ordered on the field used to interpolate.

    \author James Peachey, HEASARC
*/
#include "tip/LinearInterp.h"
#include "tip/TipException.h"

namespace tip {

  LinearInterp::LinearInterp(const Table::ConstIterator & begin, const Table::ConstIterator & end): m_begin(begin), m_end(end),
    m_match(), m_coeff() {}

  void LinearInterp::interpolate(const std::string & field, double value) {
    double x1;
    Table::ConstIterator after = m_begin;

    m_match[0] = m_end;
    m_match[1] = m_end;

    // Find first interval which might contain the value.
    for (; after != m_end; ++after) {
      (*after)[field].get(x1);
      if (x1 > value) break;
    }

    // Check ranges.
    if (m_begin == after)
      throw TipException("LinearInterp::interpolate() called for a value before the first value in range");
    else if (m_end == after)
      throw TipException("LinearInterp::interpolate() called for a value after the last value in range");

    // Success: save value of after.
    m_match[1] = after;

    // Need preceding value too.
    m_match[0] = after;
    --m_match[0];
    double x0;
    (*m_match[0])[field].get(x0);

    double D = x1 - x0;

    if (0. == D)
      throw TipException("LinearInterp::interpolate() called for an interpolation interval of 0 width");

    // Determine coefficients for interpolating all requested other fields.
    m_coeff[0] = (x1 - value) / D;
    m_coeff[1] = (value - x0) / D;

  }


  double LinearInterp::get(const std::string & field) const {
    return m_coeff[0] * (*m_match[0])[field].get() + m_coeff[1] * (*m_match[1])[field].get();
  }

  void LinearInterp::get(const std::string & field, std::vector<double> & value) const {
    std::vector<double> y0;
    std::vector<double> y1;

    // Get values from table.
    (*m_match[0])[field].get(y0);
    (*m_match[1])[field].get(y1);

    value.resize(y0.size());

    // Compute output and assign it to the output vector.
    for (std::vector<double>::size_type ii = 0; ii < y0.size(); ++ii)
      value[ii] = m_coeff[0] * y0[ii] + m_coeff[1] * y1[ii];
  }

}
