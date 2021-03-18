/** \file FileSummary.h
    \brief Class encapsulating essential information about an extension.
    \author James Peachey, HEASARC
*/
#ifndef tip_FileSummary_h
#define tip_FileSummary_h

#include <string>
#include <vector>

namespace tip {
  /** \class ExtSummary
      \brief Class encapsulating essential information about an extension.
  */
  class ExtSummary {
    public:
      /** \brief Create ExtSummary object for the given extension.
          \param ext_id The id of the extension.
      */
      ExtSummary(const std::string & ext_id);

      /** \brief Return the id of this extension.
      */
      const std::string & getExtId() const;

    private:
      std::string m_ext_id;
  };

  typedef std::vector<ExtSummary> FileSummary;

}

#endif
