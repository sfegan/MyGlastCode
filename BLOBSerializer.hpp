//-*-mode:c++; mode:font-lock;-*-

#ifndef BINARYSERIALIZE_HPP
#define BINARYSERIALIZE_HPP

#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdint.h>

#ifndef __STDC_LIMIT_MACROS
#if 0
#define __STDC_LIMIT_MACROS
#else
#error __STDC_LIMIT_MACROS and __STDC_CONSTANT_MACROS must be defined before stdint is included for the first time
#endif
#endif

#ifndef __STDC_CONSTANT_MACROS
#if 0
#define __STDC_CONSTANT_MACROS
#else
#error __STDC_LIMIT_MACROS and __STDC_CONSTANT_MACROS must be defined before stdint is included for the first time
#endif
#endif

class BLOBSerializer
{
public:
  BLOBSerializer(const std::string& filename):
    m_delete_stream(new std::ofstream(filename.c_str())), 
    m_stream(m_delete_stream) { }
  BLOBSerializer(std::ostream& stream):
    m_delete_stream(), m_stream(&stream) { }
  ~BLOBSerializer() { delete m_delete_stream; }
  bool good() { return m_stream->good(); }
  std::ostream* stream() { return m_stream; }
  template<typename T> inline bool serialize(const T& x);
private:
  std::ostream* m_delete_stream;
  std::ostream* m_stream;
};

class BLOBUnserializer
{
public:
  BLOBUnserializer(const std::string& filename):
    m_delete_stream(new std::ifstream(filename.c_str())), 
    m_stream(m_delete_stream) { }
  BLOBUnserializer(std::istream& stream):
    m_delete_stream(), m_stream(&stream) { }
  ~BLOBUnserializer() { delete m_delete_stream; }
  bool good()  { return m_stream->good(); }
  std::istream* stream() { return m_stream; }
  template<typename T> inline bool unserialize(T& x);
private:
  std::istream* m_delete_stream;
  std::istream* m_stream;
};

template<typename T> class BLOBObjectSerializer
{
public:
  static inline bool serialize(BLOBSerializer& s, const T& x)
  { return x.serializeToBlob(s); }
  static inline bool unserialize(BLOBUnserializer &s, T& x)
  { return x.unserializeFromBlob(s); }
};

template<typename T> inline bool BLOBSerializer::serialize(const T& x)
{
  return BLOBObjectSerializer<T>::serialize(*this,x);
}

template<typename T> inline bool BLOBUnserializer::unserialize(T& x)
{
  return BLOBObjectSerializer<T>::unserialize(*this,x);
}

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
//
// Specific object serializers
//
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************

#define DefineStreamSpecializedBLOBObjectSerializer(T)			\
  template<> class BLOBObjectSerializer<T>				\
  {									\
  public:								\
    static inline bool serialize(BLOBSerializer& s, const T& x)		\
    { s.stream()->write((const char*)&x,sizeof(x)); return s.good(); }	\
    static inline bool unserialize(BLOBUnserializer& s, const T& x)	\
    { s.stream()->read((char*)&x,sizeof(x)); return s.good(); }		\
  }

DefineStreamSpecializedBLOBObjectSerializer(bool);
DefineStreamSpecializedBLOBObjectSerializer(char);
DefineStreamSpecializedBLOBObjectSerializer(uint8_t);
DefineStreamSpecializedBLOBObjectSerializer(int8_t);
DefineStreamSpecializedBLOBObjectSerializer(uint16_t);
DefineStreamSpecializedBLOBObjectSerializer(int16_t);
DefineStreamSpecializedBLOBObjectSerializer(uint32_t);
DefineStreamSpecializedBLOBObjectSerializer(int32_t);
DefineStreamSpecializedBLOBObjectSerializer(uint64_t);
DefineStreamSpecializedBLOBObjectSerializer(int64_t);
DefineStreamSpecializedBLOBObjectSerializer(float);
DefineStreamSpecializedBLOBObjectSerializer(double);
DefineStreamSpecializedBLOBObjectSerializer(long double);
DefineStreamSpecializedBLOBObjectSerializer(size_t);

template<typename T>  class BLOBObjectSerializer<std::vector<T> >
{
public:
  static inline bool serialize(BLOBSerializer& s, const std::vector<T>& x)
  {
    bool ok = true;
    size_t nx(x.size());
    ok &= s.serialize(nx);
    for(unsigned ix=0;ok && ix<nx; ix++)ok &= s.serialize(x[ix]);
    return ok;
  }
  static inline bool unserialize(BLOBUnserializer& s, std::vector<T>& x)
  {
    size_t nx(0);
    if(!s.unserialize(nx))return false;
    x.resize(nx);
    for(unsigned ix=0;ix<nx; ix++)
      if(!s.unserialize(x[ix])){ x.clear(); return false; }
    return true;
  }
};

template<typename T1, typename T2> 
class BLOBObjectSerializer<std::pair<T1,T2> >
{
public:
  static inline bool serialize(BLOBSerializer& s, const std::pair<T1,T2>& x)
  { return s.serialize(x.first) && s.serialize(x.second); }
  static inline bool unserialize(BLOBUnserializer& s, std::pair<T1,T2>& x)
  { return s.unserialize(x.first) && s.unserialize(x.second); }
};

template<> class BLOBObjectSerializer<std::string>
{
public:
  static inline bool serialize(BLOBSerializer& s, const std::string& x)
  {
    size_t nx(x.length());
    if(!s.serialize(nx))return false;
    s.stream()->write(x.data(),nx);
    return s.good();
  }
  static inline bool unserialize(BLOBUnserializer& s, std::string& x)
  {
    size_t nx(0);
    if(!s.unserialize(nx))return false;
    x.clear();
    x.reserve(nx);
    for(unsigned ix=0;ix<nx; ix++)
      {
	char c;
	if(!s.unserialize(c)){ x.clear(); return false; }
	x.push_back(c);
      }
    return true;
  }
};


#endif // BINARYSERIALIZE_HPP
