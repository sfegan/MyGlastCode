/** \file IFileSvc.cxx

    \brief Factory for handling table objects.

    \author James Peachey, HEASARC
*/

#include <fstream>
#include <memory>

#include "FitsFileManager.h"
#include "FitsImage.h"
#include "FitsTable.h"
#include "FitsTipFile.h"
#include "tip/Extension.h"
#include "tip/FileSummary.h"
#include "tip/IFileSvc.h"
#include "tip/Image.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#ifndef BUILD_WITHOUT_ROOT
#include "RootTable.h"
#endif

namespace {

  std::string & s_getTmpFileName() {
    static std::string s_tmp_file_name;
    return s_tmp_file_name;
  }

}

namespace tip {

  // Deprecated.
  IFileSvc & IFileSvc::getSvc() {
    return instance();
  }

  // Get the instance of thie factory.
  IFileSvc & IFileSvc::instance() {
    static bool init_success = false;
    if (!init_success) {
      // Make sure global initialization has been performed:
      init_success = globalInit();
    }

    // Create the singleton factory.
    static IFileSvc s_file_factory;

    return s_file_factory;
  }

  // Perform global initializations.
  bool IFileSvc::globalInit() {
    bool init_success = true;
#ifndef BUILD_WITHOUT_ROOT
    init_success = RootTable::init();
#endif
    return init_success;
  }

  std::string IFileSvc::getTmpFileName() {
    return s_getTmpFileName();
  }

  void IFileSvc::setTmpFileName(const std::string & tmp_file_name) {
    s_getTmpFileName() = tmp_file_name;
  }

  // Destructor for a file service.
  IFileSvc::~IFileSvc() {
  }

  // Create a file, for now FITS only.
  void IFileSvc::createFile(const std::string & file_name, const std::string & template_name, bool clobber) {
    FitsFileManager::createFile(file_name, template_name, clobber);
  }

  TipFile IFileSvc::createMemFile(const std::string & file_name, const std::string & template_name, bool clobber) {
    return FitsFileManager::createMemFile(file_name, template_name, clobber);
  }

  TipFile IFileSvc::openFile(const std::string & file_name) {
    return new FitsTipFile(file_name);
  }

  void IFileSvc::appendImage(const std::string & file_name, const std::string & image_name,
    const ImageBase::PixelCoordinate & dims) {
    FitsFileManager::appendImage(file_name, image_name, dims);
  }

  void IFileSvc::appendTable(const std::string & file_name, const std::string & table_name) {
    FitsFileManager::appendTable(file_name, table_name);
  }

  // TODO 11: read/edit Extension/Table/Image is getting cumbersome; lots of similar methods,
  // duplicated code. In addition, error messages for files which cannot be opened tend to be
  // confusing. For instance, if bozo.fits is not present, one sees two errors: Fits file bozo.fits
  // can't be opened, and Root file bozo.fits can't be opened.
  // So: 1) Add a method which determines whether a file exists, and what type it is, and change
  // editTable so it only tries to open the correct file type. 2) Find some reasonable way to merge
  // the various edit/read methods so that there is no duplicate code, and preferably so that there
  // are fewer methods overall.

  // 5/25/2004: DONE 11: Added fileExists() method to confirm existence, and added FITS-and-Root
  // specific tests isValid(string) to FitsFileManager and RootExtensionManager, respectively.

  // Open read-write an extension in a file, be it FITS or Root, table or image.
  Extension * IFileSvc::editExtension(const std::string & file_name, const std::string & ext_name,
    const std::string & filter) {
    Extension * ext = 0;
    try {
      ext = editTable(file_name, ext_name, filter);
    } catch (const TipException & table_x) {
      try {
        ext = editImage(file_name, ext_name, filter);
      } catch (const TipException & image_x) {
        throw TipException(std::string("Could not edit extension as a table or an image:\n") + table_x.what() + "\n" +
          image_x.what());
      }
    }
    return ext;
    /*TODO 1: 4/2/2004: Memory management problem: Extension is base of Table. Extension
    has a IExtensionData and ~Extension deletes it. Currently editTable creates
    the IExtensionData and passes it to Table::Table(...) which passes it to
    Extension::Extension(...). If something throws along the way, catch 22: If
    Extension throws, editTable should delete the IExtensionData because
    ~Extension wont be called. If Table throws, editTable shouldn't delete it
    because ~Extension *will* be called. FOR NOW: take out editTable's delete,
    which may cause a memory leak in case of error, but will at least not cause
    a seg fault. */

    /* DONE 1: 4/21/2004: This is not an issue at present, because Table::Table doesn't
    throw under any circumstance. */
  }

  // Edit a image in a file, be it FITS or Root.
  TypedImage<double> * IFileSvc::editImageDbl(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    TypedImage<double> * image = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      image = new FitsTypedImage<double>(file_name, table_name, filter, false);
    else if (file_type == "root")
      throw TipException("Root images are not supported.");
    return image;
  }

  // Edit a image in a file, be it FITS or Root.
  TypedImage<float> * IFileSvc::editImageFlt(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    TypedImage<float> * image = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      image = new FitsTypedImage<float>(file_name, table_name, filter, false);
    else if (file_type == "root")
      throw TipException("Root images are not supported.");
    return image;
  }

  // Edit a image in a file, be it FITS or Root.
  TypedImage<int> * IFileSvc::editImageInt(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    TypedImage<int> * image = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      image = new FitsTypedImage<int>(file_name, table_name, filter, false);
    else if (file_type == "root")
      throw TipException("Root images are not supported.");
    return image;
  }

  // Edit a table in a file, be it FITS or Root.
  Table * IFileSvc::editTable(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    Table * table = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      table = new FitsTable(file_name, table_name, filter, false);
#ifndef BUILD_WITHOUT_ROOT
    else if (file_type == "root")
      table = new RootTable(file_name, table_name, filter, false);
#endif
    return table;
  }

  // Read-only an extension in a file, be it FITS or Root, table or image.
  const Extension * IFileSvc::readExtension(const std::string & file_name, const std::string & ext_name,
    const std::string & filter) {
    const Extension * ext = 0;
    try {
      ext = readTable(file_name, ext_name, filter);
    } catch (const TipException & table_x) {
      try {
        ext = readImage(file_name, ext_name, filter);
      } catch (const TipException & image_x) {
        throw TipException(std::string("Could not read extension as a table or an image:\n") + table_x.what() + "\n" +
          image_x.what());
      }
    }
    return ext;
  }

  // Read-only an image in a file, be it FITS or Root.
  const TypedImage<double> * IFileSvc::readImageDbl(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    TypedImage<double> * image = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      image = new FitsTypedImage<double>(file_name, table_name, filter, true);
#ifndef BUILD_WITHOUT_ROOT
    else if (file_type == "root")
      throw TipException("Root images are not supported.");
#endif
    return image;
  }

  // Read-only an image in a file, be it FITS or Root.
  const TypedImage<float> * IFileSvc::readImageFlt(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    TypedImage<float> * image = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      image = new FitsTypedImage<float>(file_name, table_name, filter, true);
#ifndef BUILD_WITHOUT_ROOT
    else if (file_type == "root")
      throw TipException("Root images are not supported.");
#endif
    return image;
  }

  // Read-only an image in a file, be it FITS or Root.
  const TypedImage<int> * IFileSvc::readImageInt(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    TypedImage<int> * image = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      image = new FitsTypedImage<int>(file_name, table_name, filter, true);
#ifndef BUILD_WITHOUT_ROOT
    else if (file_type == "root")
      throw TipException("Root images are not supported.");
#endif
    return image;
  }

  // Read-only a table in a file, be it FITS or Root.
  const Table * IFileSvc::readTable(const std::string & file_name, const std::string & table_name,
    const std::string & filter) {
    Table * table = 0;
    std::string file_type = classifyFile(file_name);
    if (file_type == "fits")
      table = new FitsTable(file_name, table_name, filter, true);
#ifndef BUILD_WITHOUT_ROOT
    else if (file_type == "root")
      table = new RootTable(file_name, table_name, filter, true);
#endif
    return table;
  }

  void IFileSvc::getFileSummary(const std::string & file_name, FileSummary & summary) {
    FitsFileManager::getFileSummary(file_name, summary);
  }

  bool IFileSvc::fileExists(const std::string & file_name) {
    std::ifstream file(file_name.c_str());
    if (!file.is_open()) return false;
    return true;
  }

  void IFileSvc::updateKeywords(const std::string & file_name, const Header::KeyValCont_t & kwds) {
    FileSummary summary;
    getFileSummary(file_name, summary);
    for (FileSummary::iterator itor = summary.begin(); itor != summary.end(); ++itor) {
      std::auto_ptr<Extension> ext(editExtension(file_name, itor->getExtId()));
      ext->getHeader().update(kwds);
    }
  }

  std::string IFileSvc::classifyFile(const std::string & file_name) {
    std::string file_type = "unknown";

    // Test whether file is a FITS file. This should be done first in case the file_name argument itself
    // has any filtering expression in it. If it does, fileExists() will say it doesn't exist even if it does.
    if (FitsFileManager::isValid(file_name)) {
      file_type = "fits";
#ifndef BUILD_WITHOUT_ROOT
    } else if (RootTable::isValid(file_name)) {
      file_type = "root";
#endif
    } else if (fileExists(file_name)) {
#ifndef BUILD_WITHOUT_ROOT
      throw TipException(std::string("File not in FITS or Root format: ") + file_name);
#else
      throw TipException(std::string("File not in FITS format: ") + file_name);
#endif
    } else {
      throw TipException(std::string("File not found: ") + file_name);
    }
    return file_type;
  }

  // Protected constructor which adds the current object to the registry of IFileSvc objects.
  IFileSvc::IFileSvc() {}

}
