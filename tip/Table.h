/** \file Table.h

    \brief High level encapsulation of tabular data.

    \author James Peachey, HEASARC
*/
#ifndef tip_Table_h
#define tip_Table_h

#include <cassert>
#include <map>
#include <string>

#include "tip/Extension.h"
#include "tip/IColumn.h"
#include "tip/Iterator.h"
#include "tip/ScalarAdaptor.h"
#include "tip/TipException.h"
#include "tip/VectorAdaptor.h"
#include "tip/tip_types.h"

namespace tip {

  class ConstTableRecord;
  class Table;

  /** \class TableCell

      \brief Encapsulation of a single table cell, which may contain (in principle) any type of object.
  */
  class TableCell {
    public:
      /** \brief Construct a TableCell object associated with a particular field in the given ConstTableRecord.
          \param record The referent ConstTableRecord object.
          \param field The name of this TableCell (the field in the ConstTableRecord).
      */
      TableCell(ConstTableRecord & record, const std::string & field): m_record(record), m_field(field),
        m_field_index(s_field_unknown) {}

      virtual ~TableCell() {}

      TableCell & operator =(const TableCell & cell);

      /** \brief Get a single value from this TableCell at the current iterator position.
          The type of the converted value is given by the template parameter.
          \param value The current value.
      */
      template <typename T>
      void get(T & value) const;

      /** \brief Get a sequence of values from this TableCell at the current iterator position.
          The type of the converted sequence is given by the template parameter.
          \param src_begin Index of the first input element.
          \param src_end Index of the element after the last input element.
          \param dest_begin Pointer to the first element in the destination sequence. The caller
          is responsible for making sure adequate space is allocated.
      */
      template <typename T>
      void get(Index_t src_begin, Index_t src_end, T * dest_begin) const;

      /** \brief Get a single value from this TableCell at the current iterator position as a double.
          (One can also use the templated get(double & value) above).
      */
      double get() const;

      /** \brief Set a sequence of values in this TableCell at the current iterator position.
          The type of the converted sequence is given by the template parameter.
          \param src_begin Pointer to the first input element.
          \param src_end Pointer to one element past the last input element.
          \param dest_begin Index of the first element in the destination sequence in the table.
      */
      template <typename T>
      void set(const T * src_begin, const T * src_end, Index_t dest_begin);

      /** \brief Set a single value in this TableCell at the current iterator position.
          The type of the source value is given by the template parameter.
          \param value The current value.
      */
      template <typename T>
      void set(const T & value);

      /** \brief Get the number of elements in the Cell.
      */
      Index_t getNumElements() const;

      /** \brief Set the number of elements in the Cell.
          \param num_elements The new number of elements.
      */
      void setNumElements(Index_t num_elements);

      /** \brief Return true if the cell contains a null (indef) value, false otherwise.
      */
      bool isNull() const;

    private:
      FieldIndex_t getFieldIndex() const;
      static const FieldIndex_t s_field_unknown = -1;
      ConstTableRecord & m_record;
      std::string m_field;
      mutable FieldIndex_t m_field_index;
  };

  /** \class ConstTableRecord

      \brief Encapsulation of a single table record.
  */
  class ConstTableRecord {
    public:
      typedef std::map<std::string, TableCell> CellCont_t;

      /** \brief Construct a ConstTableRecord object, without immediate association with a tabular data object
          Such an association may be formed later by assignment.
      */
      ConstTableRecord(): m_cells(), m_tab_data(0), m_index(0) {}

      /** \brief Construct a ConstTableRecord object, without immediate association with a tabular data object
          Such an association may be formed later by assignment.
      */
      ConstTableRecord(const ConstTableRecord & rec): m_cells(), m_tab_data(rec.m_tab_data), m_index(rec.m_index) {}

      /** \brief Construct a ConstTableRecord object, associated with the given data object and record index.
      */
      ConstTableRecord(Table * tab_data, Index_t index): m_cells(), m_tab_data(tab_data), m_index(index) {}

      /** \brief Assignment operator. Note that this behaves somewhat unusually!

          This assignment does not change this ConstTableRecord's container of TableCell objects.
          This is necessary in order to preserve the selected content of this ConstTableRecord. The
          idea is that assignment changes which data object and record index this ConstTableRecord points to
          but does not affect which fields were selected in this ConstTableRecord. This way client code
          can be certain of the continued validity of references to TableCells contained in this ConstTableRecord.
      */
    private:
      ConstTableRecord & operator =(const ConstTableRecord & rec);
    public:

      /** \brief Return a const TableCell object for the given field. The TableCell object will be created
          if it does not already exist.

          Note that there is no check at this point whether the underlying data object actually
          has a field with this name.
          \param field The name of the TableCell object (the field in this ConstTableRecord).
      */
      const TableCell & operator [](const std::string & field) const { return find_or_make(field); }

      // Client code should not normally need to call methods below here, but they are
      // public because cell and iterator abstractions need to call them.
      // Get the current record index and tab_data.
      Index_t getIndex() const { return m_index; }

      // Get the current Table pointer.
      Table * getExtensionData() { assert(m_tab_data); return m_tab_data; }

      // Get the current Table pointer.
      const Table * getExtensionData() const { assert(m_tab_data); return m_tab_data; }

      // The following methods are needed to support access to this class by the
      // iterator class.
      ConstTableRecord & itorAssign(const ConstTableRecord & record)
        { m_tab_data = record.m_tab_data; m_index = record.m_index; return *this; }
      ConstTableRecord & itorNext() { ++m_index; return *this; }
      ConstTableRecord & itorPrev() { --m_index; return *this; }
      bool itorEquals(const ConstTableRecord & record) const
        { return m_tab_data == record.m_tab_data && m_index == record.m_index; }
      bool itorLessThan(const ConstTableRecord & record) const
        { return m_tab_data == record.m_tab_data && m_index < record.m_index; }
      bool itorGreaterThan(const ConstTableRecord & record) const
        { return m_tab_data == record.m_tab_data && m_index > record.m_index; }

    protected:
      TableCell & find_or_make(const std::string & field) const;

      CellCont_t m_cells;
      Table * m_tab_data;
      Index_t m_index;
  };

  /** \class TableRecord

      \brief Encapsulation of a single table record.
  */
  class TableRecord : public ConstTableRecord {
    public:
      typedef ConstTableRecord::CellCont_t CellCont_t;

      /** \brief Construct a TableRecord object, without immediate association with a tabular data object
          Such an association may be formed later by assignment.
      */
      TableRecord(): ConstTableRecord() {}

      /** \brief Construct a TableRecord object, associated with the given tabular data object and record index.
      */
      TableRecord(Table * tab_data, Index_t index): ConstTableRecord(tab_data, index) {}

      /** \brief Copy construct a TableRecord object.
      */
      TableRecord(const TableRecord & rec): ConstTableRecord(rec) {}

      /** \brief Assignment operator. Note that this behaves somewhat unusually!

          This assignment does not change this TableRecord's container of TableCell objects.
          This is necessary in order to preserve the selected content of this TableRecord. The
          idea is that assignment changes which data object and record index this TableRecord points to
          but does not affect which fields were selected in this TableRecord. This way client code
          can be certain of the continued validity of references to TableCells contained in this TableRecord.
      */
      TableRecord & operator =(const TableRecord & rec);

      /** \brief Assignment operator. Note that this behaves somewhat unusually!

          This assignment does not change this TableRecord's container of TableCell objects.
          This is necessary in order to preserve the selected content of this TableRecord. The
          idea is that assignment changes which data object and record index this TableRecord points to
          but does not affect which fields were selected in this TableRecord. This way client code
          can be certain of the continued validity of references to TableCells contained in this TableRecord.
      */
      TableRecord & operator =(const ConstTableRecord & rec);

      /** \brief Return a TableCell object for the given field. The TableCell object will be created
          if it does not already exist.

          Note that there is no check at this point whether the underlying data object actually
          has a field with this name.
          \param field The name of the TableCell object (the field in this TableRecord).
      */
      TableCell & operator [](const std::string & field) { return find_or_make(field); }

      /** \brief Return a const TableCell object for the given field. The TableCell object will be created
          if it does not already exist.

          Note that there is no check at this point whether the underlying data object actually
          has a field with this name. Note also that this duplicates the method in the base class
          but this is necessary because the above operator [] hides the base class method.
          \param field The name of the TableCell object (the field in this ConstTableRecord).
      */
      const TableCell & operator [](const std::string & field) const { return find_or_make(field); }
  };
}

namespace tip {

  /** \class Table

      \brief High level encapsulation of tabular data, viewed as a container of Records. Access is mainly
      attained through STL-style iterators.
  */
  class Table : public Extension {
    public:
      /** \brief Helper type: encapsulation of a cell of data in this table.
      */
      typedef TableCell Cell;

      /** \brief Helper type: encapsulation of a record, containing cells of data in this table.
      */
      typedef TableRecord Record;

      /** \brief Helper type: table iterator. This currently does not work for random access.
          See TODO 7 in the source code.
      */
      typedef RandomAccessIterator<TableRecord, IndexDiff_t> Iterator;

      /** \brief Helper type: const encapsulation of a record, containing cells of data in this table.
      */
      typedef ConstTableRecord ConstRecord;

      /** \brief Helper type: const table iterator. This currently does not work for random access.
          See TODO 7 in the source code.
      */
      typedef RandomAccessIterator<TableRecord, IndexDiff_t, ConstTableRecord *, ConstTableRecord &> ConstIterator;

      /** \brief Container of field names.
      */
      typedef std::vector<std::string> FieldCont;

      /** \brief Helper type: auxilliary Cell access through a type which behaves like a primitive,
          but is connected to the table cell.
          WARNING: This class is deprecated. Don't start using it!
      */
      template <typename T>
      class Scalar : public ScalarAdaptor<T, Cell> {
        public:
          Scalar(Cell & cell): ScalarAdaptor<T, Cell>(cell) {}
          Scalar & operator =(const T & data) { ScalarAdaptor<T, Cell>::operator =(data); return *this; }
      };

      /** \brief Helper type: auxilliary Cell access through a type which behaves like a primitive,
          but is connected to the table cell.
          WARNING: This class is deprecated. Don't start using it!
      */
      template <typename T>
      class Vector : public VectorAdaptor<T, Cell> {
        public:
          Vector(Cell & cell): VectorAdaptor<T, Cell>(cell) {}
      };

      /** \brief Destructor. Closes table if it is open.
      */
      virtual ~Table() {}

      /** \brief Return an iterator pointing to the first record in the table.
      */
      Iterator begin() { return Iterator(Record(this, 0)); }

      /** \brief Return an iterator pointing past the last record in the table.
      */
      Iterator end() { return Iterator(Record(this, getNumRecords())); }

      /** \brief Return an iterator pointing to the first record in the table.
      */
// TODO Get rid of const_cast
      ConstIterator begin() const { return ConstIterator(Record(const_cast<Table *>(this), 0)); }

      /** \brief Return an iterator pointing past the last record in the table.
      */
// TODO Get rid of const_cast
      ConstIterator end() const { return ConstIterator(Record(const_cast<Table *>(this), getNumRecords())); }

      /** \brief Return the number of records contained in this table.
      */
      virtual Index_t getNumRecords() const = 0;

      /* \brief Change the number of records in the current table, adding or deleting
         rows as needed.
         \param num_records The new value for the number of records in the table.
      */
      virtual void setNumRecords(Index_t num_records) = 0;

      /** \brief Return a container of all field names valid for this table. Note that these names
          will automatically be converted to lower case, and all lookup of field names is case
          insensitive.
      */
      virtual const Table::FieldCont & getValidFields() const = 0;

      virtual IColumn * getColumn(FieldIndex_t field_index) = 0;

      virtual const IColumn * getColumn(FieldIndex_t field_index) const = 0;

      /** \brief Get an index associated with the given field (column) name.
          \param field_name The name of the field.
      */
      virtual FieldIndex_t getFieldIndex(const std::string & field_name) const = 0;

      /** \brief Copy a cell from a source extension data object to a cell in this object.
          \param src_ext The source extension data object.
          \param src_field The field identifier in the source data object.
          \param src_record The record identifier in the source data object.
          \param dest_field The field identifier in this object (the destination data object).
          \param dest_record The record identifier in this object (the destination data object).
      */
      virtual void copyCell(const Table * src_ext, FieldIndex_t src_field, Index_t src_record, FieldIndex_t dest_field,
        Index_t dest_record) = 0;

      /** \brief Copy a record from a source extension data object to a cell in this object.
          \param src_ext The source extension data object.
          \param src_record The record identifier in the source data object.
          \param dest_record The record identifier in this object (the destination data object).
      */
      virtual void copyRecord(const Table * src_ext, Index_t src_record, Index_t dest_record) = 0;

      /** \brief Append a field to the table. This will fail if a field of the same name (case insensitive) already exists.
          \param field_name The name of the field to append.
          \param format The format of the field to append, e.g. 1D for scalar double, 8J for vector long, etc.
           See Cfitsio documentation for details.
      */
      virtual void appendField(const std::string & field_name, const std::string & format) = 0;

      /** \brief Select rows in current table which match the given filtering criteria.
                 Note that this actualy changes the underlying table.
          \param filter The string containing the filtering expression.
      */
      virtual void filterRows(const std::string & filter) = 0;

  };

  /* TODO 7: 4/2/2004: 2 problems with random access: 1. operator * needs to return a
     reference to Record, but because there isn't an actual array of Records, it's not
     easy to return such a reference without leaking memory. 2. Even if 1 is solved,
     random access to Table wont be in amortized constant time, because Table may be
     seeking between calls to widely separated records (rows). Worst case: for some order
     of input, a sort could be O(N^2).
  */
}

namespace tip {
  // TableCell
  inline TableCell & TableCell::operator =(const TableCell & cell) {
    if (this != &cell) {
      const Table * src_ext = cell.m_record.getExtensionData();
      Table * dest_ext = m_record.getExtensionData();
      dest_ext->copyCell(src_ext, cell.getFieldIndex(), cell.m_record.getIndex(), getFieldIndex(), m_record.getIndex());
    }
    return *this;
  }

  // get method overloads:
  template <typename T>
  inline void TableCell::get(T & value) const {
    m_record.getExtensionData()->getColumn(getFieldIndex())->get(m_record.getIndex(), value);
  }

  template <typename T>
  inline void TableCell::get(Index_t src_begin, Index_t src_end, T * dest_begin) const {
    std::vector<T> tmp_dest(src_end);
    m_record.getExtensionData()->getColumn(getFieldIndex())->get(m_record.getIndex(), tmp_dest);
    for (IndexDiff_t ii = 0; ii < src_end - src_begin; ++ii) dest_begin[ii] = tmp_dest[src_begin + ii];
  }

  inline double TableCell::get() const {
    double retval;
    get(retval);
    return retval;
  }

  // set method overloads:
  template <typename T>
  inline void TableCell::set(const T & value) {
    m_record.getExtensionData()->getColumn(getFieldIndex())->set(m_record.getIndex(), value);
  }

  template <typename T>
  inline void TableCell::set(const T * src_begin, const T * src_end, Index_t dest_begin) {
    if (0 != dest_begin) throw TipException("BIIIIIIIIIIGGG PROBLEM");
    std::vector<T> tmp_src(src_begin, src_end);
    m_record.getExtensionData()->getColumn(getFieldIndex())->set(m_record.getIndex(), tmp_src);
  }

  inline Index_t TableCell::getNumElements() const {
    return m_record.getExtensionData()->getColumn(getFieldIndex())->getNumElements(m_record.getIndex());
  }

  inline void TableCell::setNumElements(Index_t num_elements) {
    m_record.getExtensionData()->getColumn(getFieldIndex())->setNumElements(num_elements);
  }

  inline bool TableCell::isNull() const {
    return m_record.getExtensionData()->getColumn(getFieldIndex())->isNull(m_record.getIndex());
  }

  inline FieldIndex_t TableCell::getFieldIndex() const {
    if (s_field_unknown == m_field_index)
      m_field_index = m_record.getExtensionData()->getFieldIndex(m_field);
    return m_field_index;
  }

  // ConstTableRecord
  inline TableCell & ConstTableRecord::find_or_make(const std::string & field) const {
    ConstTableRecord & rec = const_cast<ConstTableRecord &>(*this);
    CellCont_t::iterator itor = rec.m_cells.find(field);
    if (rec.m_cells.end() == itor) itor = rec.m_cells.insert(itor, std::make_pair(field, TableCell(rec, field)));
    return itor->second;
  }

  // TableRecord
  inline TableRecord & TableRecord::operator =(const TableRecord & rec) {
    if (this != &rec) {
      m_tab_data->copyRecord(rec.getExtensionData(), rec.m_index, m_index);
    }
    return *this;
  }

  inline TableRecord & TableRecord::operator =(const ConstTableRecord & rec) {
    if (this != &rec) {
      m_tab_data->copyRecord(rec.getExtensionData(), rec.getIndex(), m_index);
    }
    return *this;
  }

}

#endif
