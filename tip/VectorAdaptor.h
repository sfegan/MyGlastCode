/** \file VectorAdaptor.h

    \brief Utilities to streamline access to vector-valued cells of data inside a table.
    WARNING: This class is deprecated. Don't start using it!

    \author James Peachey, HEASARC
*/
#ifndef tip_VectorAdaptor_h
#define tip_VectorAdaptor_h

#include <cassert>

#include "tip/tip_types.h"

namespace tip {

  /** \class VectorAdaptor

      \brief Adaptor class which provides a convenient get/set interface to referents of data inside a table.
      Client code should not normally need to use this directly, but only specific subclasses of it.
      WARNING: This class is deprecated. Don't start using it!
  */
  template <typename T, typename Referent>
  class VectorAdaptor {
    public:
      class Entry {
        public:
          Entry(VectorAdaptor * adaptor, Index_t current_index): m_adaptor(adaptor), m_index(current_index) {}
          Entry & operator =(const T & data) { m_adaptor->set(m_index, data); return *this; }
          operator const T &() const { return m_adaptor->get(m_index); }
          Index_t getIndex() { return m_index; }
          void setIndex(Index_t current_index) { m_index = current_index; }
          Entry & itorNext() { ++m_index; return *this; }
          Entry & itorPrev() { --m_index; return *this; }
          bool itorEquals(const Entry & entry) const
            { return m_adaptor == entry.m_adaptor && m_index == entry.m_index; }
          bool itorLessThan(const Entry & entry) const
            { return m_adaptor == entry.m_adaptor && m_index < entry.m_index; }
          bool itorGreaterThan(const Entry & entry) const
            { return m_adaptor == entry.m_adaptor && m_index > entry.m_index; }
          void setAdaptor(VectorAdaptor * adaptor) { m_adaptor = adaptor; }

        private:
          VectorAdaptor * m_adaptor;
          Index_t m_index;
      };

//      typedef RandomAccessIterator<Entry, IndexDiff_t> Iterator;

      VectorAdaptor(Referent & referent): m_referent(&referent), m_entry(0, 0), m_begin(0), m_end(0) {
        m_entry.setAdaptor(this);
      }

      // Need to fix this:
      VectorAdaptor(const VectorAdaptor & adaptor): m_referent(adaptor.m_referent), m_entry(0, 0), m_begin(0),
        m_end(0) { m_entry.setAdaptor(this); assert(0); }

      ~VectorAdaptor() { delete [] m_begin; }

      // Need to fix this:
      VectorAdaptor & operator =(const VectorAdaptor & adaptor) { assert(0); return *this; }

      Entry & operator [](Index_t entry_index) { return getEntry(entry_index); }

      const Entry & operator [](Index_t entry_index) const
        { VectorAdaptor & self = const_cast<VectorAdaptor &>(*this); return self.getEntry(entry_index); }

      void allocate();

      const T & get(Index_t entry_index);

      Entry & getEntry(Index_t entry_index);

      Index_t getNumElements() const { return m_referent->getNumElements(); }

      void set(Index_t entry_index, const T & value);

    private:
      Referent * m_referent;
      Entry m_entry;
      T * m_begin;
      T * m_end;
  };

  template <typename T, typename Referent>
  inline void VectorAdaptor<T, Referent>::allocate() {
    Index_t vec_size = m_referent->getNumElements();
    m_begin = new T[vec_size];
    m_end = m_begin + vec_size;
  }

  template <typename T, typename Referent>
  inline const T & VectorAdaptor<T, Referent>::get(Index_t entry_index) {
    if (0 == m_begin) allocate();
    m_referent->get(0, m_end - m_begin, m_begin);
    return m_begin[entry_index];
  }

  template <typename T, typename Referent>
  inline typename VectorAdaptor<T, Referent>::Entry & VectorAdaptor<T, Referent>::getEntry(Index_t entry_index) {
    if (0 == m_begin) allocate();
    m_referent->get(0, m_end - m_begin, m_begin);
    m_entry.setIndex(entry_index);
    return m_entry;
  }

  template <typename T, typename Referent>
  inline void VectorAdaptor<T, Referent>::set(Index_t entry_index, const T & value) {
    if (0 == m_begin) allocate();
    m_begin[entry_index] = value;
    m_referent->set(m_begin, m_end, 0);
  }

}

#endif
