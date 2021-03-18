/** \file TipFile.h
    \brief Declaration for class TipFile.
    \author James Peachey, HEASARC
*/
#ifndef tip_TipFile_h
#define tip_TipFile_h

#include <string>

#include "tip/Image.h"

namespace tip {

  class Extension;
  class Table;

  /** \class ITipFile
      \brief Interface encapsulating access to an entire data file.
  */
  class ITipFile {
    public:
      /** \brief Create an ITipFile object associated with the given file.
          \param file_name The file name.
      */
      ITipFile(const std::string & file_name);

      virtual ~ITipFile() {}

      /** \brief Open the given extension (table or image) in the current file.
          \param ext_name The name of the extension to open.
      */
      virtual Extension * editExtension(const std::string & ext_name) = 0;

      /** \brief Open the given image in the current file.
          \param image_name The name of the extension to open.
      */
      virtual Image * editImage(const std::string & image_name) = 0;

      /** \brief Open the given table in the current file.
          \param table_name The name of the table to open.
      */
      virtual Table * editTable(const std::string & table_name) = 0;

      /** \brief Copy the current file to a new file.
          \param new_file_name The name of the new file to write.
          \param clobber Flag determining whether to overwrite existing files.
      */
      virtual void copyFile(const std::string & new_file_name, bool clobber = true) const = 0;

      /** \brief Clone this object.
      */
      virtual ITipFile * clone() const = 0;

      /// \brief Return the file name.
      const std::string & getName() const;

    protected:
      /** \brief Set the name of the file.
          \param file_name The new file name.
      */
      void setName(const std::string & file_name);

    private:
      std::string m_file_name;
  };

  /** \class TipFile
      \brief File format-independent data file accessor, which forwards all its methods to an ITipFile member.
  */
  class TipFile : public ITipFile {
    public:
      TipFile();

      /** \brief Copy-construct TipFile object. The underlying member is cloned.
          \param tip_file The original object being copied.
      */
      TipFile(const TipFile & tip_file);

      /** \brief Cconstruct TipFile object which forwards its members to the given ITipFile object.
                 The argument is adopted by this TipFile, and will be deleted by this TipFile's destructor.
          \param itip_file The ITipFile object to which methods will be forwarded.
      */
      TipFile(ITipFile * itip_file);

      /** \brief Destruct this TipFile and its underlying ITipFile.
      */
      virtual ~TipFile();

      /** \brief Assignment, which clones the source tip_file's underlying ITipFile.
          \param tip_file The object on the right-hand side of the assignment.
      */
      TipFile & operator =(const TipFile & tip_file);

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
      virtual void copyFile(const std::string & new_file_name, bool clobber = true) const;

      /** \brief Clone this object.
      */
      virtual ITipFile * clone() const;

    private:
      ITipFile * m_itip_file;
  };

}

#endif
