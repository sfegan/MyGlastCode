/** \file RootTable.h

    \brief Utilities to help manage Root specific table access. These classes are not part of the API.

    \author James Peachey, HEASARC
*/
#ifndef tip_RootTable_h
#define tip_RootTable_h

#include <map>
#include <string>
#include <vector>

#include "RootHeader.h"
#include "tip/IColumn.h"
#include "tip/Table.h"
#include "tip/TipException.h"
#include "tip/tip_types.h"

class TBranch;
class TFile;
class TLeaf;
class TTree;

namespace tip {

  /** \class RootTable

      \brief Low level interface to Root format extensions. This is not part of the API.

      This class is a standalone utility class which encapsulates Root access. It also
      acts as a factory for creating Root-specific header and data objects, which refer back to the
      RootTable object which created them.
  */
  class RootTable : public Table {
    public:
      /** \brief Perform global initializations needed for Root. This calls resetSigHandlers.
      */
      static bool init();

      /** \brief Reset Root's signale handlers so that Root wont interfere with debugging.
      */
      static bool resetSigHandlers();

      /** \brief Determine whether given file is a Root file.
          \param file_name The name of the file.
      */
      static bool isValid(const std::string & file_name);

      /** \brief Create an object to provide low-level access to the given Root extension.
          \param file_name The name of the Root file.
          \param ext_name The name of the Root extension.
          \param filter Root compliant filtering expression.
      */
      RootTable(const std::string & file_name, const std::string & ext_name,
        const std::string & filter = "", bool read_only = true);

      /** \brief Destructor. Closes table if it is open.
      */
      virtual ~RootTable();

      // General support for Root files:
      /** \brief Open the Root file.
      */
      void open();

      /** \brief Close the Root file.
      */
      void close();

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

      // Non-virtual helper functions for ITabularData interface:

      /** \brief Return the number of records in the current tabular data object (the number of rows
          in the Root file).
      */
      virtual Index_t getNumRecords() const;

      /* \brief Change the number of records in the current table, adding or deleting
         rows as needed.
         \param num_records The new value for the number of records in the table.
      */
      virtual void setNumRecords(Index_t num_records);

      /** \brief Return a container of all field names valid for this table:
      */
      virtual const Table::FieldCont & getValidFields() const;

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

    private:
      std::string formatWhat(const std::string & msg) const;

      RootHeader m_header;
      std::string m_file_name;
      std::string m_ext_name;
      std::string m_filter;
      std::string m_tmp_file_name;
      mutable std::map<std::string, FieldIndex_t> m_branch_lookup;
      mutable std::vector<IColumn *> m_leaves;
      Table::FieldCont m_fields;
      Index_t m_num_records;
      TFile * m_fp;
      mutable TTree * m_tree;
  };

}

#endif
