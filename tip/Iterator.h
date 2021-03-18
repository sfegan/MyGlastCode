/** \file Iterator.h

    \brief Templated iterators which iterate using a single member object of the parameterized
    type. The class pointed to by the iterator must therefore define members to implement the
    actual iterator behavior. De facto, the object is its own iterator, but it does not have
    the iterator interface, so it may define whatever operators it needs for its own purposes.

    \author James Peachey, HEASARC
*/
#ifndef tip_Iterator_h
#define tip_Iterator_h

#include <iterator>

namespace tip {

  /** \class RandomAccessIterator
      \brief Standard implementation of an iterator for situations in which the underlying container
      manages the actual data internally, and iterative access is provided through a single member object
      in the iterator, which is capable of querying the container.

      Instantiate this template to produce an iterator which may be returned by begin/end
      \param T Template argument giving the type stored in this iterator.
      \param pointer The pointer type returned by this iterator.
      \param reference The reference type returned by this iterator.
  */
  template <typename T, typename difference_t, typename pointer_t = T *, typename reference_t = T &>
  class RandomAccessIterator {
    public:
      typedef RandomAccessIterator Itor;

      // So iterator_traits will work:
      typedef std::random_access_iterator_tag iterator_category;
      typedef T value_type;
      typedef difference_t difference_type;
      typedef pointer_t pointer;
      typedef reference_t reference;

      RandomAccessIterator(): m_data() {}
      RandomAccessIterator(const T & data): m_data(data) {}

      RandomAccessIterator & operator =(const RandomAccessIterator & itor) { m_data.itorAssign(itor.m_data); return *this; }

      // In the following, the itorUpdate call is a placeholder. If e.g. Table
      // ever has a way to fill an entire row, this is how that would be accomplished.
      reference operator *() const { /* m_data.itorUpdate(); */ return m_data; }
      pointer operator ->() const { /* m_data.itorUpdate(); */ return &m_data; }

      Itor & operator ++() { m_data.itorNext(); return *this; }
      Itor operator ++(int) { Itor tmp = *this; m_data.itorNext(); return tmp; }

      Itor & operator --() { m_data.itorPrev(); return *this; }
      Itor operator --(int) { Itor tmp = *this; m_data.itorPrev(); return tmp; }

      bool operator ==(const Itor & itor) const { return m_data.itorEquals(itor.m_data); }
      bool operator !=(const Itor & itor) const { return !m_data.itorEquals(itor.m_data); }
      bool operator <(const Itor & itor) const { return m_data.itorLessThan(itor.m_data); }
      bool operator <=(const Itor & itor) const { return !m_data.itorGreaterThan(itor.m_data); }
      bool operator >(const Itor & itor) const { return m_data.itorGreaterThan(itor.m_data); }
      bool operator >=(const Itor & itor) const { return !m_data.itorLessThan(itor.m_data); }

      // The following operators pertain only to random access iterators.
      reference operator [](difference_type diff) const { return m_data.itorGet(diff); }

      Itor operator +(difference_type diff) const { Itor tmp = m_data.itorPlus(diff); return tmp; }
      Itor & operator +=(difference_type diff) { m_data.itorPlusEquals(diff); return *this; }

      Itor operator -(difference_type diff) const { Itor tmp = m_data.itorPlus(-diff); return tmp; }
      Itor & operator -=(difference_type diff) { m_data.itorPlusEquals(-diff); return *this; }

    private:
      mutable T m_data;
  };

}

#endif
