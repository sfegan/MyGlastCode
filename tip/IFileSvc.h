/** \file IFileSvc.h

    \brief Abstract factory for creating, editing and reading data objects.

    \author James Peachey, HEASARC
*/
#ifndef tip_IFileSvc_h
#define tip_IFileSvc_h

#include <string>
#include <vector>

#include "tip/FileSummary.h"
#include "tip/Header.h"
#include "tip/Image.h"
#include "tip/TipFile.h" 

namespace tip {

  class Extension;
  class Table;

  /** \class IFileSvc

      \brief Singleton factory for creating, editing and reading tables and images from files.
  */
  class IFileSvc {
    public:
      /** \brief Singleton access to I/O service objects. Deprecated: used instance() instead!
      */
      static IFileSvc & getSvc();

      /** \brief Singleton access to I/O service objects.
      */
      static IFileSvc & instance();

      /** \brief Perform initializations which are necessary at startup, mainly to handle tweaks to
          Root's global variables.
      */
      static bool globalInit();

      /// \brief Get the name of a temporary file to use when opening a large filtered file.
      static std::string getTmpFileName();

      /** \brief Set the name of a temporary file to use when opening a large filtered file.
                 It is recommended that an absolute path is used so that tip can find the temporary
                 file and delete it regardless of the current working directory.

                 Note: this currently has effect only for Root files. When keeping open multiple
                 Root tables at a time, this must be called with a distinct name before each table is opened.
          \param tmp_file_name The name of the temporary file.
      */
      static void setTmpFileName(const std::string & tmp_file_name);

      /** \brief Destruct an I/O service object.
      */
      virtual ~IFileSvc();

      /** \brief Use a FITS template to create a new file. Clobber existing files.
          \param file_name The name of the new file.
          \param template_name The name of the template file.
          \param clobber Should existing files be overwritten?
      */
      virtual void createFile(const std::string & file_name, const std::string & template_name = "", bool clobber = true);

      /** \brief Use a FITS template to create a new file in memory.
          \param file_name The name of the new file.
          \param template_name The name of the template file.
          \param clobber Should existing files be overwritten?
      */
      virtual TipFile createMemFile(const std::string & file_name, const std::string & template_name = "", bool clobber = true);

      /** \brief Open a TipFile object for access to the file as a whole.
          \param file_name The name of the file to open.
      */
      virtual TipFile openFile(const std::string & file_name);

      /** \brief Append a new image extension in a file. If the file does not exist, it will be created with
                 its primary image extension named with the image name.
          \param file_name The name of the new file.
          \param image_name The name of the new image extension.
          \param dims Set of dimensions of each axis of the image.
      */
      virtual void appendImage(const std::string & file_name, const std::string & image_name,
        const ImageBase::PixelCoordinate & dims);

      /** \brief Append a new table extension in a file. If the file does not exist, it will be created with
                 an empty primary image extension.
          \param file_name The name of the new file.
          \param table_name The name of the new table extension.
      */
      virtual void appendTable(const std::string & file_name, const std::string & table_name);

      /** \brief Open an existing extension with modification access. The actual object returned
          may be a subclass of Extension, depending on whether the object is a table or image extension.
          \param file_name The name of the file (any supported format OK).
          \param ext_name The name of the extension.
          \param filter Filtering string.
      */
      virtual Extension * editExtension(const std::string & file_name, const std::string & ext_name,
        const std::string & filter = "");

      /** \brief Open an existing image with modification access. Each pixel is treated as a float.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual Image * editImage(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "") { return editImageFlt(file_name, table_name, filter); }

      /** \brief Open an existing image with modification access. Each pixel is treated as a double.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual TypedImage<double> * editImageDbl(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing image with modification access. Each pixel is treated as a float.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual TypedImage<float> * editImageFlt(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing image with modification access. Each pixel is treated as an int.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual TypedImage<int> * editImageInt(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing table with modification access.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual Table * editTable(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing extension without modification access. The actual object returned
          may be a subclass of Extension, depending on whether the object is a table or image extension.
          \param file_name The name of the file (any supported format OK).
          \param ext_name The name of the extension.
          \param filter Filtering string.
      */
      virtual const Extension * readExtension(const std::string & file_name, const std::string & ext_name,
        const std::string & filter = "");

      /** \brief Open an existing image without modification access. Each pixel is treated as a float.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual const Image * readImage(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "") { return readImageFlt(file_name, table_name, filter); }

      /** \brief Open an existing image without modification access. Each pixel is treated as a double.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual const TypedImage<double> * readImageDbl(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing image without modification access. Each pixel is treated as a float.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual const TypedImage<float> * readImageFlt(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing image without modification access. Each pixel is treated as a double.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual const TypedImage<int> * readImageInt(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Open an existing table without modification access.
          \param file_name The name of the file (any supported format OK).
          \param table_name The name of the table.
          \param filter Filtering string.
      */
      virtual const Table * readTable(const std::string & file_name, const std::string & table_name,
        const std::string & filter = "");

      /** \brief Obtain summary of the file's contents.
          \param file_name The name of the file.
          \param summary The summary object to fill.
      */
      virtual void getFileSummary(const std::string & file_name, FileSummary & summary);

      /** \brief Check for existence of the named file.
          \param file_name The file name.
      */
      virtual bool fileExists(const std::string & file_name);

      virtual void updateKeywords(const std::string & file_name, const Header::KeyValCont_t & kwds);

    protected:
      /** \brief For singleton pattern, limit creation of IFileSvc objects to derived classes.
      */
      IFileSvc();

      std::string classifyFile(const std::string & file_name);

    private:
      // Copying file service objects is not supported.
      IFileSvc(const IFileSvc &) {}
  };

}

#endif
