/* 
  
   SourceLibrary.hpp - Stephen Fegan 
                     - sfegan@llr.in2p3.fr
                     - 26 February 2010

   Class to read Fermi source library XML file

   $Id: SourceLibrary.hpp 1952 2010-07-08 07:17:12Z sfegan $

*/

#ifndef SOURCE_LIBRARY_HPP
#define SOURCE_LIBRARY_HPP

#include<sstream>
#include<string>
#include<map>

class SourceLibrary
{
public:
  class Parameter
  {
  public:
    Parameter(const std::string& _name = "",
	      double _value = 0, double _scale = 1,
	      double _lo_bound = 0, double _hi_bound = 0, bool _free = false,
	      double _error = 0):
      m_name(_name), m_scale(_scale), 
      m_lo_bound(_lo_bound), m_hi_bound(_hi_bound),
      m_value(_value), m_free(_free), m_error(_error) { }

    const std::string& name() const { return m_name; }
    double scale() const { return m_scale; }
    double hiBound() const { return m_hi_bound; }
    double loBound() const { return m_lo_bound; }
    double value() const { return m_value; }
    double error() const { return m_error; }
    bool isFree() const { return m_free; }

    void setName(const std::string& x) { m_name = x; }
    void setScale(double x) { m_scale = x; }
    void setHiBound(double x) { m_hi_bound = x; }
    void setLoBound(double x) { m_lo_bound = x; }
    void setValue(double x) { m_value = x; }
    void setError(double x) { m_error = x; }
    void setIsFree(bool x) { m_free = x; }
    
  private:
    std::string                          m_name;
    double                               m_scale;
    double                               m_lo_bound;
    double                               m_hi_bound;
    double                               m_value;
    bool                                 m_free;
    double                               m_error;
  };

  class ParameterCollection
  {
  public:
    ParameterCollection(const std::string& _name, 
			const std::string& _type = ""):
      m_name(_name), m_type(_type), m_attributes(), m_parameters() { }

    // ACCESSER

    const std::string& name() const { return m_name; }
    const std::string& type() const { return m_type; }

    const std::map<std::string,std::string>& attributes() const 
    { return m_attributes; }
    const std::string* attribute(const std::string& name) const 
    { std::map<std::string,std::string>::const_iterator x = 
	m_attributes.find(name);
      if(x!=m_attributes.end())return &x->second; else return 0; }
    template<typename T> bool attributeAs(const std::string& name, T& x) const
    { const std::string* attr_str = attribute(name); if(attr_str == 0)return 0;
      std::istringstream(*attr_str) >> x; return true; }

    const std::map<std::string,Parameter>& parameters() const
    { return m_parameters; }
    const Parameter* parameter(const std::string& name) const 
    { std::map<std::string,Parameter>::const_iterator x = 
	m_parameters.find(name);
      if(x!=m_parameters.end())return &x->second; else return 0; }

    // SETTER

    void setType(const std::string& x) { m_type = x; }

    std::map<std::string,std::string>& attributes() { return m_attributes; }
    std::string* attribute(const std::string& name)
    { std::map<std::string,std::string>::iterator x = 
	m_attributes.find(name);
      if(x!=m_attributes.end())return &x->second; else return 0; }

    std::map<std::string,Parameter>& parameters() { return m_parameters; }
    Parameter* parameter(const std::string& name)
    { std::map<std::string,Parameter>::iterator x = m_parameters.find(name);
      if(x!=m_parameters.end())return &x->second; else return 0; }
    
    void setAttribute(const std::string& name, const std::string& value)
    { m_attributes[name] = value; }
    void clearAttributes() { m_attributes.clear(); }
    void removeAttribute(const std::string& name) 
    { std::map<std::string,std::string>::iterator x = m_attributes.find(name);
      if(x != m_attributes.end())m_attributes.erase(x); }

    void setParameter(const Parameter& x) { m_parameters[x.name()] = x; }
    void setParameter(const std::string& name,
		      double value = 0, double scale = 1,
		      double lo_bound = 0, double hi_bound = 0,
		      bool _free = false, double error = 0)
    { setParameter(Parameter(name,value,scale,lo_bound,hi_bound,_free,error));}

    void clearParameters() { m_parameters.clear(); }
    void removeParameter(const std::string& name) 
    { std::map<std::string,Parameter>::iterator x = m_parameters.find(name);
      if(x != m_parameters.end())m_parameters.erase(x); }
    
  private:
    std::string                          m_name;
    std::string                          m_type;
    std::map<std::string,std::string>    m_attributes;
    std::map<std::string,Parameter>      m_parameters;
  };

  class Source
  {
  public:
    Source(const std::string& _name = "", const std::string& _type = ""):
      m_name(_name), m_type(_type), m_attributes(), 
      m_spectrum("spectrum"), m_spatial("spatial") { }

    // ACCESSER
    
    const std::string& name() const { return m_name; }
    const std::string& type() const { return m_type; }

    const std::map<std::string,std::string>& attributes() const 
    { return m_attributes; }
    const std::string* attribute(const std::string& name) const 
    { std::map<std::string,std::string>::const_iterator x = 
	m_attributes.find(name);
      if(x!=m_attributes.end())return &x->second; else return 0; }
    template<typename T> bool attributeAs(const std::string& name, T& x) const
    { const std::string* attr_str = attribute(name); if(attr_str == 0)return 0;
      std::istringstream(*attr_str) >> x;return true; }

    const ParameterCollection& spectrum() const { return m_spectrum; }
    const ParameterCollection& spatial() const { return m_spatial; }
    
    // SETTER

    void setName(const std::string& x) { m_name = x; }
    void setType(const std::string& x) { m_type = x; }

    std::map<std::string,std::string>& attributes() { return m_attributes; }
    std::string* attribute(const std::string& name)
    { std::map<std::string,std::string>::iterator x = m_attributes.find(name);
      if(x!=m_attributes.end())return &x->second; else return 0; }

    void setAttribute(const std::string& name, const std::string& value)
    { m_attributes[name] = value; }
    void clearAttributes() { m_attributes.clear(); }
    void removeAttribute(const std::string& name) 
    { std::map<std::string,std::string>::iterator x = m_attributes.find(name);
      if(x != m_attributes.end())m_attributes.erase(x); }

    ParameterCollection& spectrum() { return m_spectrum; }
    ParameterCollection& spatial() { return m_spatial; }

  private:
    std::string                          m_name;
    std::string                          m_type;
    std::map<std::string,std::string>    m_attributes;
    ParameterCollection                  m_spectrum;
    ParameterCollection                  m_spatial;
  };

  SourceLibrary(): m_sources() { }

  // ACCESSER

  const std::map<std::string,Source>& sources() const { return m_sources; }
  const Source* source(const std::string& name) const
  { std::map<std::string,Source>::const_iterator x = m_sources.find(name);
    if(x!=m_sources.end())return &x->second; else return 0; }
  
  // SETTER
  
  std::map<std::string,Source>& sources() { return m_sources; }
  Source* source(const std::string& name)
  { std::map<std::string,Source>::iterator x = m_sources.find(name);
    if(x!=m_sources.end())return &x->second; else return 0; }

  void setSource(const Source& s) { m_sources[s.name()] = s; }

  bool readFromXML(std::string& xml_filename);
  
private:
  std::map<std::string,Source>           m_sources;
};

#endif // ifndef SOURCE_LIBRARY_HPP
