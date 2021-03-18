/* 
  
   SourceLibrary.cpp - Stephen Fegan 
                     - sfegan@llr.in2p3.fr
                     - 26 February 2010

   Class to read Fermi source library XML file

   $Id: SourceLibrary.cpp 1952 2010-07-08 07:17:12Z sfegan $

*/

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/AbstractDOMParser.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMBuilder.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMError.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMLocator.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMAttr.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/util/XMLString.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "SourceLibrary.hpp"

using namespace XERCES_CPP_NAMESPACE;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//
// Class XS - Helper for XML strings
//
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

class XS
{
public:
  XS(const std::string& x): m_x(XMLString::transcode(x.c_str())) { }
  XS(const XMLCh* const xs): m_x(XMLString::replicate(xs)) { }
  ~XS() { XMLString::release(&m_x); }
  operator const XMLCh*() const { return m_x; }
  operator std::string() const { return str(); }
  std::string str() const
  { 
    char* cstr = XMLString::transcode(m_x); std::string str(cstr);
    XMLString::release(&cstr); return str; 
  }
private:
  XMLCh* m_x;
};

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//
// Class MyErrorHandler - Handle XML pareser errors
//
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

class MyErrorHandler : public DOMErrorHandler
{
public:
  MyErrorHandler(): fSawErrors(false) { }
  ~MyErrorHandler();
  bool getSawErrors() const { return fSawErrors; }
  bool handleError(const DOMError& domError);
  void resetErrors();

private :
  MyErrorHandler(const MyErrorHandler&);
  void operator=(const MyErrorHandler&);
  bool    fSawErrors;
};

MyErrorHandler::~MyErrorHandler() { /* nothing to see here */ }

bool MyErrorHandler::handleError(const DOMError& domError)
{
  fSawErrors = true;
  if (domError.getSeverity() == DOMError::DOM_SEVERITY_WARNING)
    std::cerr << "\nWarning at file ";
  else if (domError.getSeverity() == DOMError::DOM_SEVERITY_ERROR)
    std::cerr << "\nError at file ";
  else
    std::cerr << "\nFatal Error at file ";
  
  std::cerr
    << domError.getLocation()->getURI()
    << ", line " << domError.getLocation()->getLineNumber()
    << ", char " << domError.getLocation()->getColumnNumber()
    << "\n  Message: " << domError.getMessage() << std::endl;
  
  return true;
}

void MyErrorHandler::resetErrors()
{
    fSawErrors = false;
}

bool SourceLibrary::readFromXML(std::string& xml_filename)
{
  try
    {
      XMLPlatformUtils::Initialize();
    }
  catch (const XMLException& toCatch)
    {
      std::cerr << "Error during initialization! :\n"
		<< toCatch.getMessage() << std::endl;
      return false;
    }

  // Instantiate the DOM parser.
  static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
  DOMImplementation *impl = 
    DOMImplementationRegistry::getDOMImplementation(gLS);
  DOMBuilder        *parser = 
    ((DOMImplementationLS*)impl)->
    createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
  
  parser->setFeature(XMLUni::fgDOMNamespaces, false);
  parser->setFeature(XMLUni::fgXercesSchema, false);
  parser->setFeature(XMLUni::fgXercesSchemaFullChecking, false);
  parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);

  MyErrorHandler errorHandler;
  parser->setErrorHandler(&errorHandler);

  DOMDocument *doc = 0;

  try
    {
      parser->resetDocumentPool();
      doc = parser->parseURI(xml_filename.c_str());
    }
  catch (const XMLException& toCatch)
    {
      std::cerr << "\nError during parsing: '" << xml_filename << "'\n"
		<< "Exception message is:  \n"
		<< toCatch.getMessage() << "\n" << std::endl;
      return false;
    }
  catch (const DOMException& toCatch)
    {
      const unsigned int maxChars = 2047;
      XMLCh errText[maxChars + 1];

      std::cerr << "\nDOM Error during parsing: '" << xml_filename << "'\n"
		<< "DOMException code is:  " << toCatch.code << std::endl;
      
      if(DOMImplementation::loadDOMExceptionMsg(toCatch.code,errText,maxChars))
	std::cerr << "Message is: " << errText << std::endl;
      return false;
    }
  catch (...)
    {
      std::cerr << "\nUnexpected exception during parsing: '" << xml_filename
		<< "'\n";
      return false;
    }

  if (errorHandler.getSawErrors())
    {
      std::cout << "\nErrors occurred, no output available\n" << std::endl;
      return false;
    }

  if (doc==0) 
    {
      std::cout << "No document!\n";
      return false;
    }

  DOMNode* lib = doc;
  for (lib = lib->getFirstChild(); lib != 0; lib=lib->getNextSibling())
    if(XS(lib->getNodeName()).str() == "source_library")
      break;
      
  if(lib==0)
    {
      std::cout << "No source library found!\n";
      return false;
    }
      
  for(DOMNode* src=lib->getFirstChild();src!=0;src=src->getNextSibling())
    if(XS(src->getNodeName()).str() == "source")
      {
	DOMElement* srcel = dynamic_cast<DOMElement*>(src);
	const XMLCh* srcname = srcel->getAttribute(XS("name"));
	const XMLCh* srctype = srcel->getAttribute(XS("type"));

	if(srcname==0 || srctype==0)continue;

	Source s(XS(srcname).str(),XS(srctype).str());

	DOMNamedNodeMap* attr = src->getAttributes();
	if(attr)
	  for(unsigned iattr = 0; iattr<attr->getLength(); iattr++)
	    {
	      DOMNode* n = attr->item(iattr);
	      std::string attrname = XS(n->getNodeName()).str();
	      if(attrname == "name" || attrname == "type")continue;
	      std::string attrval = XS(n->getNodeValue()).str();
	      s.setAttribute(attrname, attrval);
	    }
		
	for(DOMNode* cpt=src->getFirstChild();cpt!=0;cpt=cpt->getNextSibling())
	  {
	    std::string cptname = XS(cpt->getNodeName()).str();
	    bool is_spectrum = false;
	    if(cptname == "spectrum")is_spectrum=true;
	    else if(cptname == "spatialModel")is_spectrum=false;
	    else continue;

	    DOMElement* cptel = dynamic_cast<DOMElement*>(cpt);

	    const XMLCh* cpttype = cptel->getAttribute(XS("type"));
	    if(cpttype == 0)continue;
	    if(is_spectrum)
	      s.spectrum().setType(XS(cpttype).str());
	    else
	      s.spatial().setType(XS(cpttype).str());

	    DOMNamedNodeMap* attr = cpt->getAttributes();
	    if(attr)
	      for(unsigned iattr = 0; iattr<attr->getLength(); iattr++)
		{
		  DOMNode* n = attr->item(iattr);
		  std::string attrname = XS(n->getNodeName()).str();
		  if(attrname == "type")continue;
		  std::string attrval = XS(n->getNodeValue()).str();
		  if(is_spectrum)
		    s.spectrum().setAttribute(attrname, attrval);
		  else
		    s.spatial().setAttribute(attrname, attrval);
		}

	    for(DOMNode* par=cpt->getFirstChild();par!=0;
		par=par->getNextSibling())
	      if(XS(par->getNodeName()).str() == "parameter")
		{
		  DOMElement* parel = dynamic_cast<DOMElement*>(par);
		  const XMLCh* parname  = parel->getAttribute(XS("name"));
		  if(parname == 0)continue;

		  const XMLCh* parvalue = parel->getAttribute(XS("value"));
		  const XMLCh* parscale = parel->getAttribute(XS("scale"));
		  const XMLCh* parmax   = parel->getAttribute(XS("max"));
		  const XMLCh* parmin   = parel->getAttribute(XS("min"));
		  const XMLCh* parfree  = parel->getAttribute(XS("free"));
		  const XMLCh* parerror = parel->getAttribute(XS("error"));

		  double value   = 0;
		  double max     = 0;
		  double min     = 0;
		  double scale   = 1;
		  double is_free = false;
		  double error   = 0;
		  
		  if(parvalue)
		    std::istringstream(XS(parvalue).str()) >> value;
		  if(parmax)
		    std::istringstream(XS(parmax).str()) >> max;
		  if(parmin)
		    std::istringstream(XS(parmin).str()) >> min;
		  if(parscale)
		    std::istringstream(XS(parscale).str()) >> scale;
		  if(parfree)
		    std::istringstream(XS(parfree).str()) >> is_free;
		  if(parerror)
		    std::istringstream(XS(parerror).str()) >> error;
		  
		  Parameter p(XS(parname).str(),
			      value,scale,min,max,is_free,error);

		  if(is_spectrum)
		    s.spectrum().setParameter(p);
		  else
		    s.spatial().setParameter(p);
		}

	    setSource(s);
	  }
      }

  parser->release();
  XMLPlatformUtils::Terminate();
  return true;
}
