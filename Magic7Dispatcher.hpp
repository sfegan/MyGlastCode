//-*-mode:c++; mode:font-lock;-*-

#ifndef MAGIC7DISPATCHER_HPP
#define MAGIC7DISPATCHER_HPP

#include <stdint.h>
#include <vector>

#include "Util.hpp"
#include "FT2.hpp"

class Magic7Dispatcher
{
public:
  Magic7Dispatcher(FITSVectorVisitor<FT2>* visitor):
    m_visitor(visitor) { }
  
  FITSVectorVisitor<FT2>* setVisitor(FITSVectorVisitor<FT2>* visitor)
  { 
    FITSVectorVisitor<FT2>* ov = m_visitor;
    m_visitor = visitor;
    return ov;
  }

  unsigned dispatchVector(const std::string& filename);

private:
  FITSVectorVisitor<FT2>* m_visitor;
};

class Magic7PacketTool
{
public:
  Magic7PacketTool():
    m_pkt_file(), m_pkt() { }
  ~Magic7PacketTool();
  
  unsigned load(const std::string& filename);
  unsigned save(std::ostream& os) const;
  unsigned save(const std::string& filename) const;

  void sort();
  unsigned deleteDuplicates();

private:

  static char* readFileToBuffer(const std::string& filename, unsigned& nread);

  struct Packet
  {
    uint8_t type;
    uint32_t tsecs;
    uint32_t tfracs;
    unsigned ifile;
    unsigned ioffset;

    bool operator<(const Packet& o) const 
    {
      return((tsecs<o.tsecs)
	     ||(tsecs==o.tsecs&&(tfracs<o.tfracs
				 ||(tfracs==o.tfracs&&type>o.type))));
    }
  };

  static unsigned pktLen(uint8_t type) 
  {
    if(type==0)return 1+2*4+6*4+2*1;
    else if(type==1)return 1+2*4+4*8+3*4;
    else return 0;
  }

  std::vector<std::string>  m_pkt_file;
  std::vector<Packet> m_pkt;
};


#endif
