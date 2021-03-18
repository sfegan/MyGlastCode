/** \file KeyRecord.h
    \brief Interface for KeyRecord class.
    \authors Lawrence Brown, HEASARC/GSSC
             James Peachey, HEASARC/GSSC
*/
#ifndef tip_KeyRecord_h
#define tip_KeyRecord_h

#include <sstream>
#include <string>

namespace tip {

  /** \class KeyRecord
      \brief Encapsulation of a keyword, considered as a record with name, value and comment fields.

             See also the Keyword class. The difference between them is that Keyword obejcts
             are bound to a specific data header, whereas KeyRecord objects simply contain a
             copy of the entire keyword record stored as a single string.
  */
  class KeyRecord {
    public:
      /** \brief Construct a key record object from the given full record.
          \param record String containing the full key record, in form name = value / comment.
      */
      KeyRecord(const std::string & record = "");
  
      /** \brief Construct a key record object from each field specified separately.
          \param name The name of the keyword.
          \param value The value of the parameter.
          \param comment The comment string.
      */
      template <typename T>
      KeyRecord(const std::string & name, const T & value, const std::string & comment);

      /// \brief Returns whether the value in the record is undefined.
      bool empty() const;

      /// \brief Retrieve the record as a string.
      const std::string & get() const;

      /** \brief Assign the given input to the record.
          \param record The input record.
      */
      void set(const std::string & record);

      /// \brief Retrieve the name of the keyword record.
      std::string getName() const;

      /// \brief Retrieve the value field of the record as a string.
      std::string getValue() const;

      template <typename T>
      void getValue(T & value) const;

      /** \brief Assign the given value to the record's value field.
          \param value The new value.
      */
      void setValue(const std::string & value);

      /** \brief Assign the given value to the record's value field.
          \param value The new value.
      */
      template <typename T>
      void setValue(const T & value);

      /// \brief Retrieve the comment of the keyword record.
      std::string getComment() const;

    private:
      std::string m_record;
  };

  template <typename T>
  inline void KeyRecord::getValue(T & value) const {
    std::string str_value = getValue();
    std::stringstream ss;
    ss.precision(24);
    ss << str_value;
    ss >> value;
  }

  template <>
  inline void KeyRecord::getValue<bool>(bool & value) const {
    std::string str_value = getValue();
    if (str_value == "T") value = true;
    else value = false;
  }

  template <typename T>
  inline KeyRecord::KeyRecord(const std::string & name, const T & value, const std::string & comment): m_record() {
    // Create blank keyword, then assign a value to it.
    std::ostringstream os;
    os.width(11);
    os << std::left << name << "/ " << comment;
    m_record = os.str();
    setValue(value);
  }

  template <typename T>
  inline void KeyRecord::setValue(const T & value) {
    std::ostringstream os;
    os.precision(15);
    os << std::showpoint << value;
    setValue(os.str());
  }

  template <>
  inline void KeyRecord::setValue<bool>(const bool & value) {
    if (value) setValue("T");
    else setValue("F");
  }

}

#endif
