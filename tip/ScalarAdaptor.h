/** \file ScalarAdaptor.h

    \brief Utilities to streamline access to cells of data inside a table.
    WARNING: This class is deprecated. Don't start using it!

    \author James Peachey, HEASARC
*/
#ifndef tip_ScalarAdaptor_h
#define tip_ScalarAdaptor_h

namespace tip {

  /** \class ScalarAdaptor

      \brief Adaptor class which provides a convenient get/set interface to referents of data inside a table.
      Client code should not normally need to use this directly, but only specific subclasses of it.
      WARNING: This class is deprecated. Don't start using it!
  */
  template <typename T, typename Referent>
  class ScalarAdaptor {
    public:
      /** \brief Construct a ScalarAdaptor object which refers to the given Referent object.
          WARNING: This class is deprecated. Don't start using it!
          \param referent The Referent object.
      */
      ScalarAdaptor(Referent & referent): m_data(), m_referent(&referent) {}

      /** \brief Assignment from Referent. This changes which Referent object this ReferenceAdpator refers to.
          WARNING: This class is deprecated. Don't start using it!
          \param referent The new referent Referent object.
      */
      ScalarAdaptor & operator =(Referent & referent) { m_referent = &referent; return *this; }

      /** \brief Assignment from templated parameter type. This will write the assigned value into the
          referent to which this object refers. This does not change the Referent object this ScalarAdaptor
          refers to.
          WARNING: This class is deprecated. Don't start using it!
          \param data The source value for the assignment.
      */
      ScalarAdaptor & operator =(const T & data);

      /** \brief Retrieve the current templated parameter data value of this object.
          WARNING: This class is deprecated. Don't start using it!
      */
      operator const T & () const { m_referent->get(const_cast<T &>(m_data)); return m_data; }

    private:
      T m_data;
      Referent * m_referent;
  };

  template <typename T, typename Referent>
  ScalarAdaptor<T, Referent> & ScalarAdaptor<T, Referent>::operator =(const T & data) {
    if (m_data != data) {
      m_data = data;
      m_referent->set(m_data);
    }
    return *this;
  }

}

#endif
