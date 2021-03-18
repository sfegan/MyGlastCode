/**
    \mainpage tip package

    \author  James Peachey peachey@lheamail.gsfc.nasa.gov

    \section intro Introduction
    Provides a generic, iterator-based interface to tabular data.
    Includes file-related abstractions which are independent of the low
    level data format (e.g. FITS or ROOT).

    <hr>
    \section info Information For Clients
    This section shows how client code can and should use the tip
    classes to gain access to tabular data. The examples used to
    illustrate these concepts are all excerpted from tip's sample.cxx
    program, which can be compiled and run (provided you have an
    appropriate input file!)

    \subsection read Reading Tabular Data
    The first step is for the client to create a data access object
    representing the table:

\verbatim
#include "tip/IFileSvc.h"
#include "tip/Table.h"

  using namespace tip;
  const Table * const_table = 0;

    // Example 1: Opening a table for read-only access.
    const_table = IFileSvc::instance().readTable("day023.fits", "LAT_Event_Summary");
\endverbatim

    Table is the name of the main tip class clients need to worry about.
    The expression on the right ("IFileSvc::instance()::readTable(...)")
    calls a function which opens the given file, selects the indicated
    extension (or TTree in Root parlance) and returns a pointer to a const Table
    object which may subsequently be used to read the data. Writing is also
    supported, as described in another example below. (Some details
    for the curious: IFileSvc is a singleton abstract factory, and instance()
    returns a reference to this factory.)

    Next, one might want to read keywords from the header of the table:

\verbatim
#include "tip/Header.h"
    // Example 2:
    const Header & header = const_table->getHeader();
    double tstart;
    header["tstart"].get(tstart);
\endverbatim

    The first line creates a reference to const_table's const header object.
    The next creates a local variable to hold the value of
    the TSTART keyword. The third line causes the value of the
    TSTART keyword in the table to be copied into the local tstart
    variable. This action may also be completed in one step, but for
    clarity it is shown here in two. (Some details: Header::operator [] looks up the
    keyword in the header's container and returns a Keyword object.
    Keyword's get() method is templated and supported for all primitive
    types.)

    Another operation of interest is to read values from one or more
    columns (or TLeafs in Root parlance):

\verbatim
    // Example 3:
    // Loop over all records (rows) and extract values of ph_time column.
    for (Table::ConstIterator itor = const_table->begin(); itor != const_table->end(); ++itor) {

      // Local double variable to hold the value of the field for each row:
      double ph_time_dbl = (*itor)["ph_time"].get();

      // Do something useful with ph_time_dbl's value here ...
    }
\endverbatim

    The first clause of the for loop declares an object (itor) which
    acts as a pointer to a sequence of rows in the table. This
    object is initialized to point to the first row of the table
    by the second half of the expression ("= const_table->begin()").
    The second clause of the for loop ("itor != const_table->end()")
    causes the loop to terminate after processing the last row.
    The last clause ("++itor") causes the iterator (itor) to go on
    to the next value. Inside the loop, the iterator is dereferenced and
    then the value of the field of interest for the current record
    is obtained as a double ("(*itor)["ph_time"].get()") and stored in
    the local variable ph_time_dbl.

    In greater detail, the ConstIterator object dereferences to an object
    of a type named ConstRecord (within the Table class namespace). The ConstRecord
    type encapsulates the concept of a single record (i.e. row) in the
    table. The ConstRecord type has an operator [](const std::string &) which
    returns a reference to a const object of a type named Cell (again, within the
    Table class namespace) representing a particular cell in the table. Finally, the
    Cell type contains a get() method to return the value as a double.
    The following alternate implementation of the for loop above makes
    these details more explicit:

\verbatim
    // Example 4:
    // Completely equivalent to Example 3, just provides explicit typenames.
    // Loop over all records (rows) and extract values of ph_time column.
    for (Table::ConstIterator itor = const_table->begin(); itor != const_table->end(); ++itor) {

      // Dereference the iterator and bind it to a local reference:
      Table::ConstRecord & record = *itor;

      // Create a local reference representing the field (ph_time_cell) of
      // interest:
      const Table::Cell & ph_time_cell(record["ph_time"]);

      // Get the current value:
      double ph_time_dbl = ph_time_cell.get();

      // Do something useful with ph_time_dbl's value here ...
    }
\endverbatim

    \subsection modify Modifying Tabular Data
    The examples above demonstrate the process of opening a table in a
    file for read-only access and then iterating through it using objects
    and methods which specifically do not modify the object representing
    the Table nor the file whence it came. It is also possible to perform
    parallel actions to these which can modify the table:

\verbatim
  Table * table = 0;

    // Example 5:
    using namespace tip;
    table = IFileSvc::instance().editTable("day023.fits", "LAT_Event_Summary");

    Header & header = table->getHeader();
    double tstart;
    header["tstart"].get(tstart);

    // Modify tstart keyword:
    double offset = 86400.;
    tstart -= offset;
    header["tstart"].set(tstart);

    // Modify ph_time column:
    // Loop over all records (rows), extract and modify values of ph_time column.
    for (Table::Iterator itor = table->begin(); itor != table->end(); ++itor) {

      // Local double variable to hold the value of the field for each row:
      double ph_time_dbl = (*itor)["ph_time"].get();

      // Modify value:
      ph_time_dbl -= offset;

      // Write it back to the table:
      (*itor)["ph_time"].set(ph_time_dbl);
    }
\endverbatim

    Note that the object used to access the table is now declared
    as type Table *, not const Table *, and it is created using
    a method called editTable, instead of readTable. A reference to
    a Header object (instead of const Header) is obtained and then
    used to modify the tstart keyword. Next an iterative loop is shown
    which uses types which ultimately provide modify access to the
    table's cells.

    <hr>
    \section notes Release Notes
    \section requirements Requirements
    There shall be one or more abstractions which manage data files.
    Abstract interfaces shall be used for the API to decouple client
    code from the low level file format. Concrete implementations for FITS
    and ROOT format files shall be provided. However, because the FITS
    format is more restricted in the data structures it can sensibly
    store, the abstract interface shall conform conceptually to the
    FITS format. Restated, there shall be no methods in the abstract
    interface which provide data in a format which cannot easily be
    stored in FITS files. Essentially this means that the file interface
    will provide data in the form of table abstractions which will be
    described below. The purpose of these data file-related classes
    is to encapsulate all low level file access, including buffering,
    reading and writing data, selecting portions of data, etc.

    Abstractions shall also be designed which represent generic
    tabular data stored in a file. As with the file-related classes, abstract
    interfaces shall be used to allow client code transparently
    to access tabular data stored in a file for any supported low level
    file format. The abstract interface shall conform as closely as possible
    to Standard Template Library-style (STL-style) containers in that it
    will provide iterator-based row-oriented access to the table.
    Column-oriented access will also be supported, perhaps as a special
    case of row access in which the row contains a single field. For
    purposes of these abstractions, a table consists of a group of key
    value pairs representing FITS-style keywords, and a group of identical
    records containing one or more fields of arbitrary type and dimension.
    Access to each of these parts, that is, to the keywords (FITS header)
    and records (FITS data unit) shall function independently. Thus, if no key-value
    pairs are present (e.g. in the case of a ROOT file), the records
    in the table shall still be accessible. The fields in the records
    correspond to FITS table columns in the case of FITS files. A FITS
    image may be represented as a table with one field: the image.
    Additional image-specific abstractions may also be needed.

    For efficiency and convenience, it shall be possible to select
    subsets of data from the table-related classes. Client code will thus
    not be required to read the whole record in order to access a particular
    subset of fields.

    The concrete implementations for file-based tables (that is, concrete
    table implementations which are provided by the file-related classes,
    and which remain linked to their underlying files) shall support
    large files by providing transparent buffering as needed for efficient
    access and conservative memory usage. In this way it shall be possible
    for client code to iterate over an entire table of data, treating it
    as if it were all present in memory without explicitly managing
    this memory.

    The file-related interface shall also provide tables which
    are decoupled from the file after they are read. Clients shall be
    able to specify a given range of record (row) numbers, and receive back
    from the file interface a table implementation which was read from
    the table, but no longer relates to it. This is for two purposes.
    First, some tables may be iterated over many times, and the performance
    penalty for repeatedly reading the same information from the file
    may be unacceptable. Second, an algorithm may need to be able to modify
    a table which was read from a read-only file. In this case, it would
    be useful to have a memory-resident, modifiable table representing
    the data which was in the file.

    The file-related interface shall also support the creation of new
    files and the modification of existing files from table objects.
    This mode of access may be stream-like.

    \section plan Design Plan
    The design will be realized in a series of steps. The first step
    is to flesh out the parts of the interfaces which pertain
    to file-based tables, and implement these interfaces for the case
    of FITS files.

    \subsection table Table Class
    This is the heart of the package, and the main interface for table access.
    It contains static methods which can open or create any of the standard
    supported table types. Currently only read-only access to existing
    FITS files is supported. However it is envisioned that Root files will
    be supported, at least for reading. This class also defines Iterator,
    Record and Cell classes. The Iterator class is mostly self-explanatory.
    Dereferencing an Iterator object produces an object of type Record. A
    Record is a container of Cells. Currently it is an associative container,
    associating each Cell with its name. The Cell name corresponds to the field
    in the row, i.e. the column name in a FITS file. The Cell contains methods
    to read and write values to/from the underlying table. (Again currently,
    only reading is supported.) Just like a FITS table, a Cell can contain
    scalar data of a primitive type or character string nature, or vector data.
    Client code can request from a Cell that it read data in any of these
    forms, and if possible, the request will be satisfied.

    \subsection scalar Table::Scalar Class
    WARNING: This class is deprecated. Don't start using it!
    Templated class which serves as an adaptor for the Cell class. It may
    be templated on any primitive type or string type (currently only type
    double will work.) A Scalar<T> object binds itself to a Cell reference, and
    provides a conversion operator to type T, and assignment from type T. This
    allows client code to use a Scalar<T> object as if it were of type T, but the
    Scalar<T> object will read/write values from the table.

    \subsection image Image Class
    The Image class offers access to FITS images (the interface could
    also support Root but this is currently not implemented.)
    There are getPixel/setPixel methods which allow read/write
    access to individual pixels in an image, and also get/set
    methods to read/write entire images as a std::vector<float>.

    Clients may obtain Image objects using IFileSvc's readImage and
    editImage methods, which are similar to readTable and editTable,
    respectively.

    \section jobOptions jobOptions
    Not applicable.

    \section todo Open Issues
\verbatim
       2. 4/2/2004: Access should also be possible by extension number as well as name.
       5. 4/2/2004: In RootExtensionManager::open() there are calls to Root's
          dynamic loader, which should be eliminated if possible by a requirements pattern.
          See RootExtensionManager:cxx: TODO 5.
       7. 4/2/2004: Table's iterator has problems with random access. See Table.h: TODO 7.
       8. 4/2/2004: Use of Record within a functor used with a for_each is problematic.
          See src/test/test_tip_main.cxx: TODO 8.
          
\endverbatim

    \section done Resolved Issues
\verbatim
       1. 4/2/2004: Memory leak in IFileSvc::editTable. 4/21/2004: See IFileSvc.cxx: DONE 1.
       3. 4/2/2004: For read-only access, need const Table * IFileSvc::readTable() and
          typedef Table::ConstIterator begin() const etc. Implemented in version v1r4p0.
       4. 4/2/2004: This document's examples include some mistakes! They should be redone
          based on current live sample code. Done 4/6/2004.
       6. 4/2/2004: Add method to Table to facilitate "browsing", e.g. getting a list of
          fields in the table (as distinct from whatever happens to be in a Record).
          4/21/2004: This capability was added in v1r3p0.
       9. 4/2/2004: Bug in certain versions of cfitsio may prevent filtering syntax from
          working. See FitsExtensionManager.cxx: DONE 9.
      10. 4/21/2004: Add capability to insert a field (column) in a table. Done 4/23/2004.
      11. 5/17/2004:  Refactor IFileSvc and FitsFileManager/FitsExtensionManager to remove
          redundancy and improve error messages. See src/IFileSvc.cxx: TODO 11. Done 5/25/2004.
\endverbatim

*/
