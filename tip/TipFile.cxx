#include "tip/TipFile.h"

namespace tip {

  ITipFile::ITipFile(const std::string & name): m_file_name(name) {}

  const std::string & ITipFile::getName() const { return m_file_name; }

  void ITipFile::setName(const std::string & file_name) { m_file_name = file_name; }

  TipFile::TipFile(): ITipFile(""), m_itip_file(0) {}

  TipFile::TipFile(const TipFile & tip_file): ITipFile(tip_file.getName()), m_itip_file(tip_file.m_itip_file->clone()) {}

  TipFile::TipFile(ITipFile * itip_file): ITipFile(itip_file->getName()), m_itip_file(itip_file) {}

  TipFile::~TipFile() { delete m_itip_file; }

  TipFile & TipFile::operator =(const TipFile & tip_file) {
    delete m_itip_file;
    m_itip_file = tip_file.m_itip_file->clone();
    setName(tip_file.getName());
    return *this;
  }

  Extension * TipFile::editExtension(const std::string & ext_name) { return m_itip_file->editExtension(ext_name); }

  Image * TipFile::editImage(const std::string & image_name) { return m_itip_file->editImage(image_name); }

  Table * TipFile::editTable(const std::string & table_name) { return m_itip_file->editTable(table_name); }

  void TipFile::copyFile(const std::string & new_file_name, bool clobber) const { m_itip_file->copyFile(new_file_name, clobber); }

  ITipFile * TipFile::clone() const {
    // Pointer to be returned.
    ITipFile * itip_file = 0;

    // See if this object's ITipFile object is itself another TipFile.
    TipFile * tf = dynamic_cast<TipFile *>(m_itip_file);
    if (0 == tf) {
      // Underlying object is not a TipFile, so it is safe to clone it and return a new TipFile which owns
      // the clone.
      itip_file = new TipFile(m_itip_file->clone());
    } else {
      // Underlying object is a TipFile, so clone its underlying object rather than directly cloning it.
      itip_file = tf->m_itip_file->clone();
    }

    return itip_file;
  }

}
