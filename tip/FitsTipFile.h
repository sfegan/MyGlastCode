/** \file FitsTipFile.h
    \brief Declaration for class FitsTipFile.
    \author James Peachey, HEASARC
*/
#ifndef tip_FitsTipFile_h
#define tip_FitsTipFile_h

#include "fitsio.h"
#include "tip/TipFile.h"

namespace tip {

  /** \class FitsTipFile
      \brief File format-independent data file accessor, which forwards all its methods to an IFitsTipFile member.
  */
  class FitsTipFile : public ITipFile {
    public:
      /** \brief Create a FitsTipFile representing the given file. The file must already exist.
          \param file_name The file to be opened.
      */
      FitsTipFile(const std::string & file_name);

      /** \brief Create a FitsTipFile representing the given file. The file will be created.
          \param file_name The file to be created.
          \param template_name The template file to use.
          \param clobber If true, file will be overwritten if it exists.
      */
      FitsTipFile(const std::string & file_name, const std::string & template_name, bool clobber);

      /** \brief Copy-construct a FitsTipFile. The underlying cfitsio file will be reopened.
          \param file The source object being copied.
      */
      FitsTipFile(const FitsTipFile & file);

      virtual ~FitsTipFile();

      /** \brief Open the given extension (table or image) in the current file.
          \param ext_name The name of the extension to open.
      */
      virtual Extension * editExtension(const std::string & ext_name);

      /** \brief Open the given image in the current file.
          \param image_name The name of the extension to open.
      */
      virtual Image * editImage(const std::string & image_name);

      /** \brief Open the given table in the current file.
          \param table_name The name of the table to open.
      */
      virtual Table * editTable(const std::string & table_name);

      /** \brief Copy the current file to a new file.
          \param new_file_name The name of the new file to write.
          \param clobber Flag determining whether to overwrite existing files.
      */
      void copyFile(const std::string & new_file_name, bool clobber = true) const;

      /** \brief Clone this object.
      */
      virtual ITipFile * clone() const;

    protected:
      virtual void openFile();

      virtual void closeFile(bool update_checksum, int status);

    private:
      fitsfile * m_fp;
      bool m_read_only;
  };
}

#endif
