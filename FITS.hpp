//-*-mode:c++; mode:font-lock;-*-

#ifndef FITS_HPP
#define FITS_HPP

#include <vector>
#include <cmath>

#include "VSDataConverter.hpp"

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

// ----------------------------------------------------------------------------
// 
// FITS header
//
// ----------------------------------------------------------------------------

class FITSHeader
{
public:
  FITSHeader(): m_kv() { }
  FITSHeader(const tip::Header & header): m_kv() { loadFromHeader(header); }
  void loadFromFITS(const std::string& filename, const std::string& tablename);
  void loadFromHeader(const tip::Header & header);
  bool has(const std::string& k) const { return m_kv.find(k) != m_kv.end(); }
  unsigned numEntries(const std::string& k) const 
  {
    KV::const_iterator ikv = m_kv.find(k);
    if(ikv == m_kv.end())return 0;
    return ikv->second.size();
  }
  bool get(std::string& str, const std::string& k, unsigned ientry=0) const
  {
    KV::const_iterator ikv = m_kv.find(k);
    if(ikv == m_kv.end())return false;
    if(ientry>=ikv->second.size())return false;
    str = ikv->second[ientry];
    return true;
  }
  std::string get(const std::string& k, unsigned ientry=0)
  {
    std::string str;
    get(str,k,ientry);
    return str;
  }
  template<class T> bool getAs(T& t, const std::string& k, unsigned ientry=0)
  {
    std::string str;
    bool good = get(str, k, ientry);
    if(good)good = VERITAS::VSDataConverter::fromString(t, str);
    return good;
  }
  static std::string nk(const std::string& k, unsigned n)
  {
    std::string _k(k);
    _k += VERITAS::VSDataConverter::toString(n);
    return _k;
  }

private:
  typedef std::map<std::string, std::vector<std::string> > KV;
  KV m_kv;
};

// ----------------------------------------------------------------------------
//
// FITS visitor and dispatcher
//
// ----------------------------------------------------------------------------

template<class T> class FITSVectorVisitor
{
public:
  virtual ~FITSVectorVisitor()
  {
    // nothing to see here
  }

  virtual void visitFileVector(const std::string& filename,
			       const std::string& tablename,
			       unsigned nrow, FITSHeader& header)
  {
    // nothing to see here
  }

  virtual void leaveFileVector()
  {
    // nothing to see here
  }

  virtual void visitElement(unsigned irow, T& t)
  {
    // nothing to see here
  }
};

template<class T> class PassThroughFITSVectorVisitor:
  public FITSVectorVisitor<T>
{
public:
  PassThroughFITSVectorVisitor(FITSVectorVisitor<T>* visitor)
    : FITSVectorVisitor<T>(), m_visitor(visitor), m_nelements()
  {
    // nothing to see here
  }

  virtual ~PassThroughFITSVectorVisitor()
  {
    // nothing to see here
  }

  virtual void visitFileVector(const std::string& filename,
			       const std::string& tablename,
			       unsigned nrow, FITSHeader& header)
  {
    m_visitor->visitFileVector(filename, tablename, nrow, header);
  }

  virtual void leaveFileVector()
  {
    m_visitor->leaveFileVector();
  }

  virtual void visitElement(unsigned irow, T& t)
  {
    m_visitor->visitElement(irow, t);
    m_nelements++;
  }

  unsigned nElementsDispatched() const { return m_nelements; }

protected:
  FITSVectorVisitor<T>* m_visitor;
  unsigned m_nelements;
};

template<class T> class MultiFITSVectorVisitor:
  public FITSVectorVisitor<T>
{
public:
  MultiFITSVectorVisitor(bool copy_element = true)
    : FITSVectorVisitor<T>(), m_visitors(), m_copy_element(copy_element)
  {
    // nothing to see here
  }

  virtual ~MultiFITSVectorVisitor()
  {
    // nothing to see here
  }

  void addVisitor(FITSVectorVisitor<T>* v) { m_visitors.push_back(v); }

  virtual void visitFileVector(const std::string& filename,
			       const std::string& tablename,
			       unsigned nrow, FITSHeader& header)
  {
    for(typename std::vector<FITSVectorVisitor<T> *>::iterator ivisitor =
	  m_visitors.begin(); ivisitor != m_visitors.end(); ivisitor++)
      (*ivisitor)->visitFileVector(filename, tablename, nrow, header);
  }

  virtual void leaveFileVector()
  {
    for(typename std::vector<FITSVectorVisitor<T> *>::iterator ivisitor =
	  m_visitors.begin(); ivisitor != m_visitors.end(); ivisitor++)
      (*ivisitor)->leaveFileVector();
  }

  virtual void visitElement(unsigned irow, T& t)
  {
    for(typename std::vector<FITSVectorVisitor<T> *>::iterator ivisitor =
	  m_visitors.begin(); ivisitor != m_visitors.end(); ivisitor++)
      if(m_copy_element)
	{
	  T tt(t);
	  (*ivisitor)->visitElement(irow, tt);
	}
      else
	{
	  (*ivisitor)->visitElement(irow, t);
	}
  }

protected:
  std::vector<FITSVectorVisitor<T> *> m_visitors;
  bool m_copy_element;
};

template<class T> class DumpToOStreamVisitor:
  public FITSVectorVisitor<T>
{
public:
  DumpToOStreamVisitor(std::ostream& stream)
    : FITSVectorVisitor<T>(), m_stream(stream)
  {
    // nothing to see here
  }

  virtual ~DumpToOStreamVisitor()
  {
    // nothing to see here
  }

  virtual void visitElement(unsigned irow, T& t)
  {
    t.dumpToOStream(m_stream);
  }

protected:
  std::ostream& m_stream;
};

template<class T> class FITSVectorDispatcher
{
public:
  FITSVectorDispatcher(FITSVectorVisitor<T>* visitor):
    m_visitor(visitor) { }

  FITSVectorVisitor<T>* setVisitor(FITSVectorVisitor<T>* visitor)
  { 
    FITSVectorVisitor<T>* ov = m_visitor;
    m_visitor = visitor;
    return ov;
  }

  unsigned dispatchVector(const std::string& filename,
			  const std::string& tablename,
			  const typename T::FITSFillOptions& opt =
			  typename T::FITSFillOptions())
  {
    const tip::Table* fits =
      tip::IFileSvc::instance().readTable(filename, tablename);

    const tip::Header & fits_header(fits->getHeader());
    FITSHeader header(fits_header);
    unsigned nrow(0);
    if(header.has("NAXIS2"))header.getAs(nrow,"NAXIS2");

    m_visitor->visitFileVector(filename,tablename,nrow,header);

    tip::Table::ConstIterator it = fits->begin();
    tip::Table::ConstRecord & datum = *it;

    unsigned ntotal(0);
    for ( ; it != fits->end(); ++it, ntotal++) 
      {
	T t;
	try
	  {
	    t.fillFromFITS(datum, header, opt);
	  }
	catch(tip::TipException& x)
	  {
	    std::cerr << "TipException at record: " << ntotal << '\n';
	    throw;
	  }
	m_visitor->visitElement(ntotal, t);
      }

    m_visitor->leaveFileVector();

    return ntotal;
  }

private:
  FITSVectorVisitor<T>* m_visitor;
};

// ----------------------------------------------------------------------------
//
// Simple function to fill a vector from FITS file
//
// ----------------------------------------------------------------------------

template<class T> class FITSFillVectorVisitor: public FITSVectorVisitor<T>
{
public:
  FITSFillVectorVisitor(std::vector<T>& v): m_v(v) { }
  virtual void visitFileVector(const std::string& filename,
			       const std::string& tablename,
			       unsigned nrow){ m_v.reserve(m_v.size()+nrow); }
  virtual void visitElement(unsigned irow, T& t) { m_v.push_back(t); }
private:
  std::vector<T>& m_v;
};

template<class T> unsigned
fillVectorFromFITS(std::vector<T>& t_vec,
		   const std::string& filename, const std::string& tablename,
		   const typename T::FITSFillOptions& opt = 
		   typename T::FITSFillOptions())
{
  t_vec.clear();
  FITSFillVectorVisitor<T> visitor(t_vec);
  FITSVectorDispatcher<T> dispatcher(&visitor);
  unsigned ntotal = dispatcher.dispatchVector(filename, tablename, opt);
  return ntotal;
};

#endif
