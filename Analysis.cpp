//-*-mode:c++; mode:font-lock;-*-

#include <iostream>
#include <fstream>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "VSALinearLeastSquares.hpp"

#include "Analysis.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;
using namespace VERITAS::VSAMath;

// ----------------------------------------------------------------------------
//
// SOURCE FINDER
//
// ----------------------------------------------------------------------------

FoundSourceVisitor::~FoundSourceVisitor()
{
  // nothing to see here
}

void FoundSourceVisitor::
visitFSEvent(unsigned irow, FT1& t, 
	     unsigned isrc, double ra, double dec, double dist)
{
  // nothing to see here
}

SourceFinder::
SourceFinder(FoundSourceVisitor* visitor, double d_max,
	     const std::vector<double>& cat_ra,
	     const std::vector<double>& cat_dec):
  PassThroughFITSVectorVisitor<FT1>(visitor),
  m_visitor(visitor), m_d_max(d_max),
  m_ncat(), m_cat_ra(cat_ra), m_cat_dec(cat_dec),
  m_cat_sra(), m_cat_cra(), m_cat_sdec(), m_cat_cdec(),
  m_nfsevents()
{
  m_ncat = m_cat_ra.size();
  m_cat_sra.resize(m_ncat);
  m_cat_cra.resize(m_ncat);
  m_cat_sdec.resize(m_ncat);
  m_cat_cdec.resize(m_ncat);

  for(unsigned icat=0;icat<m_ncat;icat++)
    {
      m_cat_sra[icat] = std::sin(m_cat_ra[icat]);
      m_cat_cra[icat] = std::cos(m_cat_ra[icat]);
      m_cat_sdec[icat] = std::sin(m_cat_dec[icat]);
      m_cat_cdec[icat] = std::cos(m_cat_dec[icat]);
    }  
}

SourceFinder::~SourceFinder()
{
  // nothing to see here
}

void SourceFinder::visitElement(unsigned irow, FT1& event)
{
  PassThroughFITSVectorVisitor<FT1>::visitElement(irow, event);

  double ra = event.ra;
  double dec = event.dec;

  double sra = std::sin(ra);
  double cra = std::cos(ra);
  double sdec = std::sin(dec);
  double cdec = std::cos(dec);

  double dmin = 0;
  unsigned dmin_icat = 0;
  for(unsigned icat=0;icat<m_ncat;icat++)
    {
      double d = 
	approx_sphere_dist_nt(sra,cra,sdec,cdec,
			      m_cat_sra[icat],m_cat_cra[icat],
			      m_cat_sdec[icat],m_cat_cdec[icat]);
      if(icat==0 || d<dmin)dmin=d,dmin_icat=icat;
    }

  double cat_ra = m_cat_ra[dmin_icat];
  double cat_dec = m_cat_dec[dmin_icat];

  dmin = sphere_dist(ra, dec, cat_ra, cat_dec);
  if(dmin<m_d_max)
    {
      m_visitor->visitFSEvent(irow, event, dmin_icat, cat_ra, cat_dec, dmin);
      m_nfsevents++;
    }
}

// ----------------------------------------------------------------------------
//
// RADEC ORIENTATION
//
// ----------------------------------------------------------------------------

RADECRotation::~RADECRotation()
{
  // nothing to see here
}

void RADECRotation::visitElement(unsigned irow, FT1& event)
{
  FT1 mod_event = event;
  Vec3D d = Vec3D::makePolar(M_PI_2-mod_event.dec,mod_event.ra);
  d.rotate(m_radec_rot);
  mod_event.ra = d.phi();
  if(mod_event.ra<0)mod_event.ra += 2.0*M_PI;
  mod_event.dec = M_PI_2-d.theta();
  PassThroughFITSVectorVisitor<FT1>::visitElement(irow, mod_event);
}

// ----------------------------------------------------------------------------
//
// SPACECRAFT ORIENTATION
//
// ----------------------------------------------------------------------------

SOEventVisitor::~SOEventVisitor()
{
  // nothing to see here
}

void SOEventVisitor::visitSOEvent(unsigned irow, FT1& event, Vec3D& so)
{
  // nothing to see here
}


EventRADECRecalc::~EventRADECRecalc()
{
  // nothing to see here
}

void EventRADECRecalc::visitSOEvent(unsigned irow, FT1& event, Vec3D& so)
{
  FT1 mod_event = event;
  unsigned nfc = m_fourier_dtheta.size();
  double dtheta = 0;
  for(unsigned ifc=0;ifc<nfc;ifc++)
    dtheta += m_fourier_dtheta[ifc]*std::sin(mod_event.theta*double(ifc+1));
  mod_event.theta += dtheta;
  Vec3D d = Vec3D::makePolar(mod_event.theta,mod_event.phi);
  d.rotate(m_sc_rot);
  d.rotate(so);
  mod_event.ra = d.phi();
  mod_event.dec = M_PI_2-d.theta();
  m_visitor->visitElement(irow, mod_event);
}
 
void EventRADECRecalc::visitFileVector(const std::string& filename,
					  const std::string& tablename,
				       unsigned nrow, FITSHeader& header)
{
  m_visitor->visitFileVector(filename, tablename, nrow, header);
}

void EventRADECRecalc::leaveFileVector()
{
  m_visitor->leaveFileVector();
}

EODatEntryMatcher::
EODatEntryMatcher(SOEventVisitor* visitor, const std::string& eodat_fn)
  : PassThroughFITSVectorVisitor<FT1>(visitor), 
    m_visitor(visitor), m_eo_map(), m_nmatched()
{
  int fd = open(eodat_fn.c_str(), O_RDONLY);
  if(fd<0)throw std::string("Could not open "+eodat_fn);

  unsigned magic(0);
  if(read(fd,&magic,sizeof(magic))==sizeof(magic) && magic==0xFEEDBEEF)
    {
      unsigned runno;
      unsigned eventno;
      double x;
      double y;
      double z;
      while(read(fd,&runno,sizeof(runno))==sizeof(runno)
	    && read(fd,&eventno,sizeof(eventno))==sizeof(eventno)
	    && read(fd,&x,sizeof(x))==sizeof(x)
	    && read(fd,&y,sizeof(y))==sizeof(y)
	    && read(fd,&z,sizeof(z))==sizeof(z))
	{
	  std::pair<unsigned,unsigned> key(runno,eventno);
	  Vec3D so(x,y,z);
	  m_eo_map[key]=so;
	}
      close(fd);
    }
  else
    {
      close(fd);

      std::ifstream streamf(eodat_fn.c_str());
      std::string line;
  
      std::getline(streamf,line);
      while(streamf)
	{
	  std::istringstream streaml(line);
	  unsigned runno;
	  unsigned eventno;
	  double x;
	  double y;
	  double z;
	  streaml >> runno >> eventno >> x >> y >> z;
	  std::pair<unsigned,unsigned> key(runno,eventno);
	  Vec3D so(x,y,z);
	  m_eo_map[key]=so;
	  std::getline(streamf,line);
	}
    }
}

void EODatEntryMatcher::visitElement(unsigned irow, FT1& event)
{
  std::pair<unsigned,unsigned> key(event.run_id,event.event_id);
  std::map<std::pair<unsigned,unsigned>, VERITAS::VSAAlgebra::Vec3D>::iterator
    ifind = m_eo_map.find(key);
  if(ifind != m_eo_map.end())
    {
      m_visitor->visitSOEvent(irow, event, ifind->second);
      m_nmatched++;
    }
}

QuadInterpSOCalculator::~QuadInterpSOCalculator()
{
  // nothing to see here
}

static inline double I(double x,
		       const std::vector<double>& xi,
		       const std::vector<double>& yi)
{
  unsigned ndata(xi.size());
  Data data(ndata);
  for(unsigned idata=0;idata<ndata;idata++)
    {
      data[idata].x = xi[idata];
      data[idata].y = yi[idata];
      data[idata].sigma = 0.001;
      //std::cerr << xi[idata] << ' ' << yi[idata] << '\n';
    }

  VSAAlgebra::VecND p = PolyFit::fit(2, data);
  double y = PolyFit::val(p, x);
  //std::cerr << x << ' ' << y << '\n';
  //std::cout << p.size() << '\n';
  return y;
}
		       
void QuadInterpSOCalculator::visitElement(unsigned irow, FT1& event)
{
  PassThroughFITSVectorVisitor<FT1>::visitElement(irow, event);

  while(m_ift2!=m_ft2.end() && m_ift2->t_stop<event.time)m_ift2++;
  while(m_ift2->t_start>event.time && m_ift2!=m_ft2.begin())m_ift2--;
   
  if(m_ift2!=m_ft2.end() 
     && m_ift2->t_stop>=event.time && m_ift2->t_start<=event.time)
    {
      std::vector<FT2>::const_iterator ift2 = m_ift2;

#if 1
      double t = event.time;

      std::vector<double> ti;
      std::vector<double> scz_xi;
      std::vector<double> scz_yi;
      std::vector<double> scz_zi;
      std::vector<double> scx_xi;
      std::vector<double> scx_yi;
      std::vector<double> scx_zi;

      const double T = 10.2;
      while(ift2!=m_ft2.begin() && (t-ift2->t_start)<T)ift2--;
      while(ift2!=m_ft2.end() && (ift2->t_start-t)<T)
	{
	  Vec3D sczi = Vec3D::makePolar(M_PI_2-ift2->scz_dec, ift2->scz_ra);
	  Vec3D scxi = Vec3D::makePolar(M_PI_2-ift2->scx_dec, ift2->scx_ra);
	  ti.push_back(ift2->t_start-t);
	  scz_xi.push_back(sczi.x());
	  scz_yi.push_back(sczi.y());
	  scz_zi.push_back(sczi.z());
	  scx_xi.push_back(scxi.x());
	  scx_yi.push_back(scxi.y());
	  scx_zi.push_back(scxi.z());
	  ift2++;
	}

      Vec3D scz(I(0.0,ti,scz_xi),I(0.0,ti,scz_yi),I(0.0,ti,scz_zi));
      Vec3D scx(I(0.0,ti,scx_xi),I(0.0,ti,scx_yi),I(0.0,ti,scx_zi));
      
      Vec3D so  = Vec3D::makeRotation(scz,scx,SS_31);
      m_visitor->visitSOEvent(irow, event, so);
      m_nsoevents++;
#else
      if((ift2+1)!=m_ft2.end() &&
	 fabs(event.time-(ift2+1)->t_start)<fabs(event.time-ift2->t_start))
	ift2++;

      Vec3D scz = Vec3D::makePolar(M_PI_2-ift2->scz_dec, ift2->scz_ra);
      Vec3D scx = Vec3D::makePolar(M_PI_2-ift2->scx_dec, ift2->scx_ra);
      Vec3D so  = Vec3D::makeRotation(scz,scx,SS_31);
      m_visitor->visitSOEvent(irow, event, so);
      m_nsoevents++;
#endif
    }
}

