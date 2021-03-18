/** \file FitsImage.h

    \brief Utilities to help manage FITS specific table access. These classes are not part of the API.

    \author James Peachey, HEASARC
*/
#ifndef tip_FitsImage_h
#define tip_FitsImage_h

#include <string>
#include <vector>

#include "fitsio.h"

#include "FitsHeader.h"
#include "FitsPrimProps.h"
#include "tip/Image.h"
#include "tip/tip_types.h"

namespace tip {

  /** \class FitsTypedImage

      \brief Low level interface to FITS format extensions. This is not part of the API.

      This class is a standalone utility class which encapsulates Cfitsio (fitsfile *) access. It also
      acts as a factory for creating FITS-specific header and data objects, which refer back to the
      FitsTypedImage object which created them.
  */
  template <typename T>
  class FitsTypedImage : public TypedImage<T> {
    public:
      /** \brief Create an object to provide low-level access to the given FITS extension.
          \param file_name The name of the FITS file.
          \param ext_name The name of the FITS extension.
      */
      FitsTypedImage(const std::string & file_name, const std::string & ext_name,
        const std::string & filter = "", bool read_only = true);

      /** \brief Destructor. Closes image if it is open.
      */
      ~FitsTypedImage();

      /** \brief Retrieve Header object, which is a container of FITS-like keywords, non-const version.
      */
      Header & getHeader();

      /** \brief Retrieve Header object, which is a container of FITS-like keywords, const version.
      */
      const Header & getHeader() const;

      /** \brief Returns true if the extension is a image, false otherwise.
      */
      bool isImage() const { return true; }

      /** \brief Returns true if the extension is a table, false otherwise.
      */
      bool isTable() const { return false; }

      /** \brief Return name of this extension.
      */
      const std::string & getName() const;

      /** \brief Set name of this extension.
      */
      void setName(const std::string & name);

      /** \brief Get the dimensionality of an image.
      */
      const ImageBase::PixelCoordinate & getImageDimensions() const;

      /** \brief Set the dimensionality of an image.
          \param dims Array holding the sizes in each dimension.
      */
      void setImageDimensions(const ImageBase::PixelCoordinate & dims);

      /** \brief Get a specific pixel from an image extension.
          \param coord The coordinates of the pixel.
          \param pixel The pixel value.
      */
      void getPixel(const ImageBase::PixelCoordinate & coord, double & pixel) const;

      /** \brief Set a specific pixel in an image extension.
          \param coord The coordinates of the pixel.
          \param pixel The pixel value.
      */
      void setPixel(const ImageBase::PixelCoordinate & x, const double & pixel);

      /** \brief Return a specific pixel from an image extension.
          \param coord The coordinates of the pixel.
      */
      virtual T get(const ImageBase::PixelCoordinate & coord) const;

      /** \brief Set a specific pixel in an image extension.
          \param coord The coordinates of the pixel.
          \param pixel The pixel value.
      */
      virtual void set(const ImageBase::PixelCoordinate & coord, T pixel);

      /** \brief Get an entire image, regardless of its dimensionality, as a one-dimensional array.
          \param image The array in which to store the image.
      */
      void get(std::vector<T> & image) const;

      /** \brief Get a slice of an image as a one-dimensional array.
          \param range A container of intervals which give the range of pixels in each image dimension.
          \param image The array in which to store the image.
      */
      void get(const ImageBase::PixelCoordRange & range, std::vector<T> & image) const;

      /** \brief Get an entire image, regardless of its dimensionality, as a one-dimensional array.
          \param image The array which stores the image to be written.
      */
      void set(const std::vector<T> & image);

      /** \brief Set a slice of an image as a one-dimensional array.
          \param range A container of intervals which give the range of pixels in each image dimension.
          \param image The array which stores the slice of the image to be written.
      */
      void set(const ImageBase::PixelCoordRange & range, const std::vector<T> & image);

    protected:
      void close(int status = 0);

      fitsfile * getFp() const { return m_header.getFp(); }

      bool readOnly() const { return m_header.readOnly(); }

      /** \brief Open a FITS image. Exceptions will be thrown if the extension does not exist, or if
          the extension is not an image. Normally this is called by open()
      */
      void openImage();

    private:
      std::string formatWhat(const std::string & msg) const;

      FitsHeader m_header;
      std::string m_file_name;
      std::string m_filter;
      ImageBase::PixelCoordinate m_image_dimensions;
  };

  typedef FitsTypedImage<float> FitsImage;

  template <typename T>
  inline FitsTypedImage<T>::FitsTypedImage(const std::string & file_name, const std::string & ext_name,
    const std::string & filter, bool read_only): m_header(file_name, ext_name, filter, read_only),
    m_file_name(file_name), m_filter(filter), m_image_dimensions() { openImage(); }

  // Close file automatically while destructing.
  template <typename T>
  inline FitsTypedImage<T>::~FitsTypedImage() { close(); }

  // Close file.
  template <typename T>
  inline void FitsTypedImage<T>::close(int status) {
    m_image_dimensions.clear();
    m_header.close(status);
  }

  template <typename T>
  inline Header & FitsTypedImage<T>::getHeader() { return m_header; }

  template <typename T>
  inline const Header & FitsTypedImage<T>::getHeader() const { return m_header; }

  template <typename T>
  inline const std::string & FitsTypedImage<T>::getName() const { return m_header.getName(); }

  template <typename T>
  inline void FitsTypedImage<T>::setName(const std::string & name) { m_header.setName(name); }

  template <typename T>
  inline const ImageBase::PixelCoordinate & FitsTypedImage<T>::getImageDimensions() const {
    return m_image_dimensions;
  }

  template <typename T>
  inline void FitsTypedImage<T>::setImageDimensions(const ImageBase::PixelCoordinate & dims) {
    // Make C primitive copy of array to pass to Cfitsio.
    ImageBase::PixelCoordinate::size_type naxis = dims.size();
    ImageBase::PixelCoordinate naxes(dims.begin(), dims.end());

    int status = 0;
    int bitpix = 0;

    // Get current image type, which will not be changed by the resize operation.
    fits_get_img_type(m_header.getFp(), &bitpix, &status);
    if (0 != status) throw TipException(status, formatWhat("setImageDimensions cannot determine image type"));

    // Resize the image.
    fits_resize_img(m_header.getFp(), bitpix, naxis, &*naxes.begin(), &status);
    if (0 != status) throw TipException(status, formatWhat("setImageDimensions cannot change image dimensions"));

    // Save the dimensions in the dimension member.
    m_image_dimensions = dims;
  }

  template <typename T>
  T FitsTypedImage<T>::get(const ImageBase::PixelCoordinate & coord) const {
    int status = 0;
    // Make a copy of coordinates for cfitsio to use.
    ImageBase::PixelCoordinate cf_coord(coord.size());

    // Cfitsio starts numbering at 1 not 0.
    for (ImageBase::PixelCoordinate::size_type index = 0; index != cf_coord.size(); ++index) cf_coord[index] = coord[index] + 1;

    T array[2] = { 0, 0 };

    // Read the given pixel:
    fits_read_pix(m_header.getFp(), FitsPrimProps<T>::dataTypeCode(), &*cf_coord.begin(), 1, 0, array, 0, &status);
    if (0 != status) throw TipException(status, formatWhat("get(coordinate) could not read pixel"));

    return array[0];
  }

  template <typename T>
  void FitsTypedImage<T>::set(const ImageBase::PixelCoordinate & coord, T pixel) {
    if (m_header.readOnly()) throw TipException(formatWhat("set(coordinate, pixel) called for read-only image"));
    int status = 0;
    // Make a copy of coordinates for cfitsio to use.
    ImageBase::PixelCoordinate cf_coord(coord.size());

    // Cfitsio starts numbering at 1 not 0.
    for (ImageBase::PixelCoordinate::size_type index = 0; index != cf_coord.size(); ++index) cf_coord[index] = coord[index] + 1;

    // Copy pixel into temporary array:
    T array[2] = { pixel, 0 };

    // Write the copy to the output file:
    fits_write_pix(m_header.getFp(), FitsPrimProps<T>::dataTypeCode(), &*cf_coord.begin(), 1, array, &status);
    if (0 != status) throw TipException(status, formatWhat("set(coordinate, pixel) could not write a pixel"));
  }

  template <typename T>
  inline void FitsTypedImage<T>::getPixel(const ImageBase::PixelCoordinate & coord, double & pixel) const {
    int status = 0;
    // Make a copy of coordinates for cfitsio to use.
    ImageBase::PixelCoordinate cf_coord(coord.size());

    // Cfitsio starts numbering at 1 not 0.
    for (ImageBase::PixelCoordinate::size_type index = 0; index != cf_coord.size(); ++index) cf_coord[index] = coord[index] + 1;

    double array[2] = { 0., 0. };

    // Read the given pixel:
    fits_read_pix(m_header.getFp(), TDOUBLE, &*cf_coord.begin(), 1, 0, array, 0, &status);
    if (0 != status) throw TipException(status, formatWhat("getPixel could not read pixel as a double"));

    // Copy the value just read:
    pixel = *array;
  }

  template <typename T>
  inline void FitsTypedImage<T>::setPixel(const ImageBase::PixelCoordinate & coord, const double & pixel) {
    if (m_header.readOnly()) throw TipException(formatWhat("setPixel called for read-only image"));
    int status = 0;
    // Make a copy of coordinates for cfitsio to use.
    ImageBase::PixelCoordinate cf_coord(coord.size());

    // Cfitsio starts numbering at 1 not 0.
    for (ImageBase::PixelCoordinate::size_type index = 0; index != cf_coord.size(); ++index) cf_coord[index] = coord[index] + 1;

    // Copy pixel into temporary array:
    double array[2] = { pixel, 0. };

    // Write the copy to the output file:
    fits_write_pix(m_header.getFp(), TDOUBLE, &*cf_coord.begin(), 1, array, &status);
    if (0 != status) throw TipException(status, formatWhat("setPixel could not write a double to a pixel"));
  }

  template <typename T>
  inline void FitsTypedImage<T>::openImage() {
    // Check whether the file pointer is pointing at a table:
    if (m_header.isTable()) {
      close();
      throw TipException(formatWhat("HDU is not an image"));
    }

    int naxis = 0;
    int status = 0;

    // Get number of axes:
    fits_get_img_dim(m_header.getFp(), &naxis, &status);
    if (0 != status) throw TipException(status, formatWhat("Cannot get number of dimensions of image"));

    m_image_dimensions.clear();

    // Get naxes:
    ImageBase::PixelCoordinate naxes(naxis);
    fits_get_img_size(m_header.getFp(), naxis, &*naxes.begin(), &status);
    if (0 != status) {
      throw TipException(status, formatWhat("Cannot get dimensions of each degree of freedom of image"));
    }

    m_image_dimensions.assign(naxes.begin(), naxes.end());
  }

  template <typename T>
  inline std::string FitsTypedImage<T>::formatWhat(const std::string & msg) const {
    std::ostringstream msg_str;
    msg_str << msg;
    const std::string & ext_name(getName());
    if (!ext_name.empty()) msg_str << " in extension \"" << ext_name << '"';
    msg_str << " in file \"" << m_file_name << '"';
    return msg_str.str();
  }

  template <typename T>
  inline void FitsTypedImage<T>::get(std::vector<T> & image) const {
    int status = 0;

    // Compute overall size of image.
    PixOrd_t image_size = 1;
    for (ImageBase::PixelCoordinate::const_iterator itor = m_image_dimensions.begin(); itor != m_image_dimensions.end(); ++itor)
      image_size *= *itor;

    // Make sure image is large enough to accomodate the image.
    image.resize(image_size);

    // Starting coordinate is the first pixel in each dimension.
    ImageBase::PixelCoordinate coord(m_image_dimensions.size(), 1);

    // Get the image itself.
    fits_read_pix(m_header.getFp(), FitsPrimProps<T>::dataTypeCode(), &*coord.begin(), image_size, 0, &*image.begin(), 0, &status);
  }

  template <typename T>
  inline void FitsTypedImage<T>::get(const ImageBase::PixelCoordRange & range, std::vector<T> & image) const {
    int status = 0;

    // Create arrays which contain first and last pixel in cfitsio's indexing scheme.
    std::vector<long> fpixel(range.size());
    std::vector<long> lpixel(range.size());

    // Interpret range to get arrays of first and last pixels, as well as the total image size.
    // Note: Correct for fact that lpixel is indexed starting with 1 not 0:
    //   fpixel = begin_pixel + 1 (offset for indexing)
    // However, for lpixel there is a second correction because end_pixel is defined as one past the last pixel. So:
    //   lpixel = end_pixel + 1 (offset for indexing) - 1 (end_pixel is one pixel past last pixel) = end_pixel
    // Thus, these corrections offset, and so there is no correction for lpixel.
    PixOrd_t slice_size = 1;
    for (ImageBase::PixelCoordRange::size_type index = 0; index != range.size(); ++index) {
      fpixel[index] = range[index].first + 1; // Cfitsio indexes image coordinates starting with 1 not 0.
      lpixel[index] = range[index].second; // DO NOT add 1, because range already is 1 past the last pixel.
      slice_size *= range[index].second - range[index].first;
    }

    // Resize image array.
    image.resize(slice_size);

    // Set up the down sampling array: skip no pixels.
    std::vector<long> inc(fpixel.size(), 1);

    // Get the image.
    fits_read_subset(m_header.getFp(), FitsPrimProps<T>::dataTypeCode(), &*fpixel.begin(), &*lpixel.begin(), &*inc.begin(),
      0, &*image.begin(), 0, &status);
    if (0 != status) throw TipException(status, formatWhat("could not read image subset"));
  }

  template <typename T>
  inline void FitsTypedImage<T>::set(const std::vector<T> & image) {
    int status = 0;

    // Compute overall size of image.
    PixOrd_t image_size = 1;
    for (ImageBase::PixelCoordinate::const_iterator itor = m_image_dimensions.begin(); itor != m_image_dimensions.end(); ++itor)
      image_size *= *itor;

    // Make sure no more than the image_size elements are written.
    image_size = (image_size < PixOrd_t(image.size())) ? image_size : image.size();

    // Starting coordinate is the first pixel in each dimension.
    ImageBase::PixelCoordinate coord(m_image_dimensions.size(), 1);

    // Write the image itself.
    fits_write_pix(m_header.getFp(), FitsPrimProps<T>::dataTypeCode(), &*coord.begin(), image_size,
      const_cast<T *>(&*image.begin()), &status);
  }

  template <typename T>
  inline void FitsTypedImage<T>::set(const ImageBase::PixelCoordRange & range, const std::vector<T> & image) {
    int status = 0;

    // Create arrays which contain first and last pixel in cfitsio's indexing scheme.
    std::vector<long> fpixel(range.size());
    std::vector<long> lpixel(range.size());

    // Interpret range to get arrays of first and last pixels, as well as the total image size.
    // Note: Correct for fact that lpixel is indexed starting with 1 not 0:
    //   fpixel = begin_pixel + 1 (offset for indexing)
    // However, for lpixel there is a second correction because end_pixel is defined as one past the last pixel. So:
    //   lpixel = end_pixel + 1 (offset for indexing) - 1 (end_pixel is one pixel past last pixel) = end_pixel
    // Thus, these corrections offset, and so there is no correction for lpixel.
    for (ImageBase::PixelCoordRange::size_type index = 0; index != range.size(); ++index) {
      fpixel[index] = range[index].first + 1; // Cfitsio indexes image coordinates starting with 1 not 0.
      lpixel[index] = range[index].second; // DO NOT add 1, because range already is 1 past the last pixel.
    }

    // Write the image itself.
    fits_write_subset(m_header.getFp(), FitsPrimProps<T>::dataTypeCode(), &*fpixel.begin(), &*lpixel.begin(),
      const_cast<T *>(&*image.begin()), &status);
  }

}

#endif
