/** \file FitsFileManager.h

    \brief Utilities to help manage FITS files. These classes are not part of the API.

    \author James Peachey, HEASARC
*/
#ifndef tip_FitsFileManager_h
#define tip_FitsFileManager_h

#include <string>
#include <vector>

#include "fitsio.h"
#include "tip/FileSummary.h"
#include "tip/Image.h"
#include "tip/TipFile.h"

namespace tip {
  /** \class FitsFileManager
      \brief Low level interface to FITS files. This is not part of the API. This class is concerned
      with whole files, as distinct from other classes which handle individual file extensions.
  */
  class FitsFileManager {
    public:
      /** \brief Create a new file, using an optional FITS template file. Clobber existing files.
          \param file_name The name of the new file.
          \param template_name The name of the template file. If omitted, the file created will
                 be created with an empty primary image.
          \param clobber Should existing files be overwritten?
      */
      static void createFile(const std::string & file_name, const std::string & template_name = "", bool clobber = true);

      /** \brief Use a FITS template to create a new file in memory.
          \param file_name The name of the new file.
          \param template_name The name of the template file.
          \param clobber Should existing files be overwritten?
      */
      static TipFile createMemFile(const std::string & file_name, const std::string & template_name = "", bool clobber = true);

      /** \brief Append an image extension to a fits file.
          \param file_name The name of the file to which to append.
          \param image_name The name of the new table.
          \param dims The set of sizes for each dimension of the image.
      */
      static void appendImage(const std::string & file_name, const std::string & image_name,
        const ImageBase::PixelCoordinate & dims);

      /** \brief Append a table extension to a file.
          \param file_name The name of the file to which to append.
          \param table_name The name of the new table.
      */
      static void appendTable(const std::string & file_name, const std::string & table_name);

      /** \brief Obtain summary of contents of the given file.
          \param file_name The name of the file.
          \param summary The object holding the summary.
      */
      static void getFileSummary(const std::string & file_name, FileSummary & summary);

      /** \brief Determine if the given file name is the name of a FITS file.
          \param file_name The name of the file.
      */
      static bool isValid(const std::string & file_name);

    private:
      static fitsfile * createFile(const std::string & file_name, const std::string & image_name,
        const ImageBase::PixelCoordinate & dims);
      static fitsfile * createImage(fitsfile * fp, const std::string & file_name, const std::string & image_name,
        const ImageBase::PixelCoordinate & dims);

      // Get the extsnsion identifier (name or number).
      static void getExtId(fitsfile * fp, std::string & ext_id);

      static void closeFile(fitsfile *fp, bool update_checksum, int status);
  };

}

#endif
