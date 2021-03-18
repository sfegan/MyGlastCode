#include <vector>
#include <fstream>
#include <cmath>
#include <stdint.h>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>

#ifdef __APPLE__
#include <libkern/OSByteOrder.h>
#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)
#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)
#endif

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "VERITAS/VSAAlgebra.hpp"

#include "Magic7Dispatcher.hpp"

using namespace VERITAS::VSAAlgebra;

static uint8_t readUInt8(std::istream& s)
{
  uint8_t u8;
  s.read(reinterpret_cast<char*>(&u8),sizeof(u8));
  return u8;
}

static uint32_t readUInt32(std::istream& s)
{
  uint32_t u32;
  s.read(reinterpret_cast<char*>(&u32),sizeof(u32));
  return be32toh(u32); 
}

static uint64_t readUInt64(std::istream& s)
{
  uint64_t u64;
  s.read(reinterpret_cast<char*>(&u64),sizeof(u64));
  return be64toh(u64); 
}

static float readFlt(std::istream& s)
{
  uint32_t u32 = readUInt32(s);
  float f;
  assert(sizeof(f) == sizeof(u32));
  f = *(float*)(void*)(&u32);
  return f;
}

static double readDbl(std::istream& s)
{
  uint64_t u64 = readUInt64(s);
  double d;
  assert(sizeof(d) == sizeof(u64));
  d = *(double*)(void*)(&u64);
  return d;
}

unsigned Magic7Dispatcher::dispatchVector(const std::string& filename)
{
  std::ifstream fs(filename.c_str());
  if(!fs)
    {
      return 0;
    }
  
  FITSHeader header;
  m_visitor->visitFileVector(filename,"N/A",0,header);

  uint32_t att_tsecs(0);
  uint32_t att_tfracs(0);
  double att_qx(0);
  double att_qy(0);
  double att_qz(0);
  double att_qw(0);
  float att_wx(0);
  float att_wy(0);
  float att_wz(0);
  
  unsigned ntotal(0);
  unsigned char pt = readUInt8(fs);
  while(fs)
    {
      if(pt == 1)
	{
	  att_tsecs        = readUInt32(fs);
	  att_tfracs       = readUInt32(fs);
	  att_qx           = readDbl(fs);
	  att_qy           = readDbl(fs);
	  att_qz           = readDbl(fs);
	  att_qw           = readDbl(fs);
	  att_wx           = readFlt(fs);
	  att_wy           = readFlt(fs);
	  att_wz           = readFlt(fs);

#if 0
	  std::cerr
	    << "ATT: " << att_tsecs << ' ' << att_tfracs << ' '
	    << att_qx << ' ' << att_qy << ' ' << att_qz << ' ' << att_qw << ' '
	    << att_wx << ' ' << att_wy << ' ' << att_wz << '\n';
#endif
	}
      else if(pt == 0)
	{
	  uint32_t tsecs   = readUInt32(fs);
	  uint32_t tfracs  = readUInt32(fs);
	  float x          = readFlt(fs);
	  float y          = readFlt(fs);
	  float z          = readFlt(fs);
	  float vx         = readFlt(fs);
	  float vy         = readFlt(fs);
	  float vz         = readFlt(fs);
	  uint8_t mode     = readUInt8(fs);
	  uint8_t in_saa   = readUInt8(fs);

#if 0
	  std::cerr 
	    << "ORB: " << tsecs << ' ' << tfracs << ' '
	    << x << ' ' << y << ' ' << z << ' '
	    << vx << ' ' << vy << ' ' << vz << ' ' 
	    << mode << ' ' << in_saa << '\n';
#endif
  
	  if(tsecs==att_tsecs && tfracs==att_tfracs)
	    {
	      Vec3D r(x,y,z);
	      Vec3D v(vx,vy,vz);

	      Vec3D pole(r^v);

	      RotationVec3D qrot = 
	       RotationVec3D::makeRotationQuat(att_qw, att_qx, att_qy, att_qz);

	      Vec3D scz(0.0,0.0,1.0);
	      scz.rotate(qrot);

	      Vec3D scx(1.0,0.0,0.0);
	      scx.rotate(qrot);

	      FT2 ft2;
	    
	      ft2.t_start  = double(tsecs) + double(tfracs)*1e-6;
	      ft2.t_stop   = double(tsecs) + double(tfracs)*1e-6 + 1.0;

	      ft2.zn_ra    = r.lonPos();
	      ft2.zn_dec   = r.lat();

	      ft2.in_saa   = in_saa;

	      ft2.scz_ra   = scz.lonPos();
	      ft2.scz_dec  = scz.lat();
	      ft2.scx_ra   = scx.lonPos();
	      ft2.scx_dec  = scx.lat();
	    
	      ft2.pole_ra  = pole.lonPos();
	      ft2.pole_dec = pole.lat();

	      ft2.qso_w    = att_qw; 
	      ft2.qso_x    = att_qx;
	      ft2.qso_y    = att_qy;
	      ft2.qso_z    = att_qz;

	      m_visitor->visitElement(ntotal, ft2);
	      ntotal++;
	    }
	  else
	    std::cerr << "NO ATT entry for ORB: " 
		      << tsecs << ' ' << tfracs << ' '
		      << x << ' ' << y << ' ' << z << ' '
		      << vx << ' ' << vy << ' ' << vz << ' ' 
		      << mode << ' ' << in_saa << '\n';
	}
      assert(fs);
      pt = readUInt8(fs);
    }

  return ntotal;
}

Magic7PacketTool::~Magic7PacketTool()
{
  // nothing to see here
}

char* Magic7PacketTool::
readFileToBuffer(const std::string& filename, unsigned& nread)
{
  nread = 0;

  struct stat bstat;
  int istat = stat(filename.c_str(), &bstat);
  if(istat<0)return 0;

  nread = bstat.st_size;
  char* pb = new char[nread];

  std::ifstream fs(filename.c_str());
  if(!fs)
    {
      delete[] pb;
      nread = 0;
      return 0;
    }

  unsigned iread = 0;
  while(fs && iread<nread)
    {
      fs.read(pb+iread, std::max(nread-iread, 1024U*1024U));
      iread += fs.gcount();
    }

  if(iread != nread)
    {
      delete[] pb;
      nread = 0;
      return 0;
    }

  return pb;
}

unsigned Magic7PacketTool::load(const std::string& filename)
{

  unsigned nread = 0;
  char* pb = readFileToBuffer(filename, nread);
  if(pb == 0)return 0;
  m_pkt_file.push_back(filename);

  unsigned ifile = m_pkt_file.size()-1;
  unsigned iread = 0;
  unsigned npkt = 0;
  while(iread<nread)
    {
      Packet p;
      char* pdata = pb+iread;
      p.ifile = ifile;
      p.ioffset = iread;
      p.type = *pdata;
      iread += pktLen(p.type);
      assert(iread<=nread);
      p.tsecs = be32toh(*reinterpret_cast<uint32_t*>(pdata+1));
      p.tfracs = be32toh(*reinterpret_cast<uint32_t*>(pdata+5));
      m_pkt.push_back(p);
      npkt++;
    }

  delete[] pb;
  return npkt;
}

unsigned Magic7PacketTool::save(std::ostream& os) const
{
  unsigned npkt=m_pkt.size();
  unsigned ipkt=0;

  unsigned nread = 0;
  unsigned ifile = m_pkt_file.size();
  char* pb(0);

  while(os && ipkt<npkt)
    {
      if(m_pkt[ipkt].ifile != ifile)
	{
	  delete[] pb;
	  pb = 0;
	  ifile = m_pkt[ipkt].ifile;
	  pb = readFileToBuffer(m_pkt_file[ifile], nread);
	  assert(pb);
	}

      char* pdata = pb+m_pkt[ipkt].ioffset;
      unsigned plen = pktLen(m_pkt[ipkt].type);
      assert(m_pkt[ipkt].ioffset + plen <= nread);
      os.write(pdata, plen);
      ipkt++;
    }
  return ipkt;
}

unsigned Magic7PacketTool::save(const std::string& filename) const
{
  std::ofstream fs(filename.c_str());
  if(!fs)return 0;
  return save(fs);
}

void Magic7PacketTool::sort()
{
  std::sort(m_pkt.begin(), m_pkt.end());
}

unsigned Magic7PacketTool::deleteDuplicates()
{
  std::vector<unsigned> last_packet;
  unsigned npkt=m_pkt.size();
  unsigned ipkt=0;
  unsigned jpkt=0;

  while(ipkt<npkt)
    {
      unsigned itype = m_pkt[ipkt].type;
      if(itype>=last_packet.size())last_packet.resize(itype+1,npkt);

      if(jpkt != ipkt)m_pkt[jpkt] = m_pkt[ipkt];

      if(jpkt < last_packet[itype]
	 || m_pkt[jpkt].tsecs != m_pkt[last_packet[itype]].tsecs 
	 || m_pkt[jpkt].tfracs != m_pkt[last_packet[itype]].tfracs)
	{
	  last_packet[itype]=jpkt;
	  jpkt++;
	}

      ipkt++;
    }
  m_pkt.resize(jpkt);
  return npkt-jpkt;
}
