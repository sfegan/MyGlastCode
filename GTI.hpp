//-*-mode:c++; mode:font-lock;-*-

#ifndef GTI_HPP
#define GTI_HPP

#include <iostream>
#include <algorithm>

#include "BLOBSerializer.hpp"
#include "FITS.hpp"

class GTI
{
public:
  GTI(): t_start(), t_stop() { }
  GTI(const double t): t_start(t), t_stop(t) { }
  GTI(const double start, const double stop): t_start(start), t_stop(stop) { }

  double t_start;
  double t_stop;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  static const char* tableName() { return "gti"; }
  static unsigned loadFromFITS(std::vector<GTI>& gti_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());

  bool operator<(const GTI& o) const { return t_start<o.t_start; }

  bool serializeToBlob(BLOBSerializer& s) const
  { return s.serialize(t_start) && s.serialize(t_stop); }
  bool unserializeFromBlob(BLOBUnserializer& s)
  { return s.unserialize(t_start) && s.unserialize(t_stop); }
};

class GTIRange
{
public:
  GTIRange(): m_gti() { }

  // ITERATOR
  std::vector<GTI>::size_type size() const { return m_gti.size(); }
  bool empty() const { return m_gti.empty(); }
  typedef std::vector<GTI>::const_iterator const_iterator;
  const_iterator begin() const { return m_gti.begin(); }
  const_iterator end() const { return m_gti.end(); }
  
  // GETTERS
  unsigned nGTIs() const { return m_gti.size(); }
  const std::vector<GTI>& getGTIs() const { return m_gti; }
  double totalTime() const;

  // SETTERS
  void loadGTIsFromFITS(const std::string& filename, 
			const std::string& tablename = GTI::tableName(),
			const GTI::FITSFillOptions& opt = 
			GTI::FITSFillOptions());
  template<class T> void addGTIs(const std::vector<T>& x)
  {
    m_gti.reserve(m_gti.size() + x.size());
    for(typename std::vector<T>::const_iterator ix=x.begin();ix!=x.end();ix++)
      m_gti.push_back(GTI(ix->t_start, ix->t_stop));
    normalize();
  }
  template<class T> void setGTIs(const std::vector<T>& x)
  {
    m_gti.clear();
    addGTIs(x);
  }
  template<class T> void addGTI(const T& x)
  {
    if(empty())
      {
	m_gti.push_back(GTI(x.t_start, x.t_stop));
	m_last = begin();
      }
    else if(x.t_start>=m_gti.back().t_start)
      {
	if(x.t_start>m_gti.back().t_stop)
	  {
	    m_gti.push_back(GTI(x.t_start, x.t_stop));
	    m_last = begin();
	  }
	else if(x.t_stop>m_gti.back().t_stop)
	  m_gti.back().t_stop = x.t_stop;
      }
    else
      {
	std::cout << "addGTI - SLOW\n";
	// SLOW VERSION
	m_gti.push_back(GTI(x.t_start, x.t_stop));
	normalize();
      }
  }

  // FIND TIME
  const_iterator findGTI(double t) const
  {
    if(empty())return end();
    else if(m_last==end())
      {
	if(t>=m_gti.back().t_stop)return m_last;
      }
    else if(t>=m_last->t_start)
      {
	if(t<m_last->t_stop)return m_last;
	m_last++;
	if(m_last==end())return m_last;
	else if(t<m_last->t_start)
	  {
	    m_last--;
	    return end();
	  }
	else if(t<m_last->t_stop)return m_last;
      }
       
    m_last = std::upper_bound(begin(),end(),GTI(t));

    if(m_last==begin())return end();
    m_last--;
    if(t>=m_last->t_start && t<m_last->t_stop)return m_last;
    else return end();    
  }

  // FIND TIME
  const_iterator _findGTI(double t) const
  {
    const_iterator i1 = findGTI(t);

    const_iterator i2 = begin();
    while(i2 != end())
      {
	if((t>=i2->t_start)&&(t<i2->t_stop))break;
	i2++;
      }
    assert(i1==i2);
    return i1;
  }
  
  bool serializeToBlob(BLOBSerializer& s) const;
  bool unserializeFromBlob(BLOBUnserializer& s);

private:
  void normalize();
  std::vector<GTI> m_gti;
  mutable std::vector<GTI>::const_iterator m_last;
};

#endif
