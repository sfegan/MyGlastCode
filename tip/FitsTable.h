/** \file FitsTable.h

    \brief Utilities to help manage FITS specific table access. These classes are not part of the API.

    \author James Peachey, HEASARC
*/
#ifndef tip_FitsTable_h
#define tip_FitsTable_h

#include <map>
#include <string>
#include <vector>

#include "fitsio.h"

#include "FitsHeader.h"
#include "tip/IColumn.h"
#include "tip/Table.h"
#include "tip/tip_types.h"

namespace tip {

  /** \class FitsTable

      \brief Low level interface to FITS format extensions. This is not part of the API.

      This class is a standalone utility class which encapsulates Cfitsio (fitsfile *) access. It also
      acts as a factory for creating FITS-specific header and data objects, which refer back to the
      FitsTable object which created them.
  */
  class FitsTable : public Table {
    public:
      /** \brief Create an object to provide low-level access to the given FITS extension.
          \param file_name The name of the FITS file.
          \param ext_name The name of the FITS extension.
      */
      FitsTable(const std::string & file_name, const std::string & ext_name,
        const std::string & filter = "", bool read_only = true);

      /** \brief Destructor. Closes table if it is open.
      */
      virtual ~FitsTable();

      /** \brief Close the FITS file.
          \param status The status to use when closing the file. Note this is not a reference!
      */
      void close(int status = 0);

      /** \brief Retrieve Header object, which is a container of FITS-like keywords, non-const version.
      */
      virtual Header & getHeader();

      /** \brief Retrieve Header object, which is a container of FITS-like keywords, const version.
      */
      virtual const Header & getHeader() const;

      /** \brief Returns true if the extension is a image, false otherwise.
      */
      virtual bool isImage() const { return false; }

      /** \brief Returns true if the extension is a table, false otherwise.
      */
      virtual bool isTable() const { return true; }

      /** \brief Return name of this extension.
      */
      virtual const std::string & getName() const;

      /** \brief Set name of this extension.
      */
      virtual void setName(const std::string & name);

      /** \brief Return the number of records in the current tabular data object (the number of rows
          in the FITS file).
      */
      virtual Index_t getNumRecords() const;

      /* \brief Change the number of records in the current table, adding or deleting
         rows as needed.
         \param num_records The new value for the number of records in the table.
      */
      virtual void setNumRecords(Index_t num_records);

      /** \brief Return a container of all field names valid for this table:
      */
      virtual const FieldCont & getValidFields() const;

      virtual IColumn * getColumn(FieldIndex_t field_index);

      virtual const IColumn * getColumn(FieldIndex_t field_index) const;

      /** \brief Get an index associated with the given field (column) name.
          \param field_name The name of the field.
      */
      virtual FieldIndex_t getFieldIndex(const std::string & field_name) const;

      /** \brief Copy a cell from a source extension data object to a cell in this object.
          \param src_ext The source extension data object.
          \param src_field The field identifier in the source data object.
          \param src_record The record identifier in the source data object.
          \param dest_field The field identifier in this object (the destination data object).
          \param dest_record The record identifier in this object (the destination data object).
      */
      virtual void copyCell(const Table * src_ext, FieldIndex_t src_field, Index_t src_record, FieldIndex_t dest_field,
        Index_t dest_record);

      /** \brief Copy a record from a source extension data object to a cell in this object.
          \param src_ext The source extension data object.
          \param src_record The record identifier in the source data object.
          \param dest_record The record identifier in this object (the destination data object).
      */
      virtual void copyRecord(const Table * src_ext, Index_t src_record, Index_t dest_record);

      /** \brief Append a field to the table.
          \param field_name The name of the field to append.
          \param format The format of the field to append, e.g. 1D for scalar double, 8J for vector long, etc.
           See Cfitsio documentation for details.
      */
      virtual void appendField(const std::string & field_name, const std::string & format);

      /** \brief Select rows in current table which match the given filtering criteria.
                 Note that this actualy changes the underlying table.
          \param filter The string containing the filtering expression.
      */
      virtual void filterRows(const std::string & filter);

      fitsfile * getFp() const { return m_header.getFp(); }

      bool readOnly() const { return m_header.readOnly(); }

    protected:
      /** \brief Open the FITS table. Exceptions will be thrown if the extension does not exist, or if
          the extension is not a table. Normally this is called by open()
      */
      void openTable();

      /** \brief Get all pertinent information about a column, and store it:
          \param col_name The name of the column.
          \param col_num The number of the column in the FITS file.
      */
      void getColumnInfo(const std::string & col_name, Index_t col_num);

    private:
      std::string formatWhat(const std::string & msg) const;

      FitsHeader m_header;
      std::string m_file_name;
      std::string m_filter;
      std::map<std::string, FieldIndex_t> m_col_name_lookup;
      FieldCont m_fields;
      std::vector<IColumn *> m_columns;
      Index_t m_num_records;
  };

  // Copying cells.
  inline void FitsTable::copyCell(const Table * src_ext, FieldIndex_t src_field, Index_t src_record,
    FieldIndex_t dest_field, Index_t dest_record) {
    getColumn(dest_field)->copy(src_ext->getColumn(src_field), src_record, dest_record);
  }

  // Copying records.
  inline void FitsTable::copyRecord(const Table * src_ext, Index_t src_record, Index_t dest_record) {
    for (Table::FieldCont::iterator itor = m_fields.begin(); itor != m_fields.end(); ++itor) {
      copyCell(src_ext, src_ext->getFieldIndex(*itor), src_record, getFieldIndex(*itor), dest_record);
    }
  }
}

#endif
