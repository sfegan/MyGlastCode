//-*-mode:c++; mode:font-lock;-*-

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"

#include "FT1.hpp"
#include "FT2.hpp"
#include "Catalog.hpp"
#include "Util.hpp"
#include "MyMinuit.hpp"
#include "PSF.hpp"
#include "Analysis.hpp"
#include "Magic7Dispatcher.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

class OrbitCalculator: public FITSVectorVisitor<FT2>
{
public:
  OrbitCalculator(double t0, double tslice, unsigned nslice);
  virtual ~OrbitCalculator();
  virtual void visitElement(unsigned irow, FT2& t);
private:
  unsigned iel(unsigned ix, unsigned iy) { return iy*m_nx+ix; }
  bool writeSlice(unsigned new_islice);
  static inline double ea(double ct);

  Vec3D m_g2j;

  double m_t0;
  double m_tslice;
  unsigned m_nslice;
  unsigned m_islice;
  unsigned m_nx;
  unsigned m_ny;
  std::vector<double> m_data;
  Vec3D m_j2o;
  std::vector<double> m_x;
  std::vector<double> m_y;
};

inline double OrbitCalculator::ea(double ct)
{
  // Sum of front and back EA @ 1GeV
  static const double A[] = 
    { 0.0076896, 0.0171072, 0.0294624, 0.0427680, 0.0675648, 0.0825120,
      0.1083456, 0.1327968, 0.1559520, 0.1855872, 0.2077056, 0.2373408,
      0.2590272, 0.2858112, 0.3065472, 0.3381696, 0.3564864, 0.3874176,
      0.4078080, 0.4367520, 0.4447008, 0.4726080, 0.4892832, 0.5116608,
      0.5307552, 0.5404320, 0.5689440, 0.5971104, 0.5852736, 0.6209568,
      0.6411744, 0.7088256, 0.7088256 };

  ct -= 0.2;
  if(ct<0)return 0;
  const double q = floor(ct*40.0);
  const double x = ct-q*0.025;
  const unsigned i = q;
  return A[i]*(1.0-x)+A[i+1]*x;
}

OrbitCalculator::OrbitCalculator(double t0, double tslice, unsigned nslice)
  : FITSVectorVisitor<FT2>(), 
    m_g2j(Vec3D::makeRotationGalToJ2000()), 
    m_t0(t0), m_tslice(tslice), m_nslice(nslice), m_islice(nslice),
    m_nx(200), m_ny(100), m_data(), m_j2o(), m_x(), m_y()
{ 
  m_x.resize(m_nx);
  for(unsigned ix=0;ix<m_nx;ix++)
    m_x[ix] = ((double(ix)+0.5)/double(m_nx)*2.0-1.0)*M_PI;
  m_y.resize(m_ny);
  for(unsigned iy=0;iy<m_ny;iy++)
    m_y[iy] = std::asin((double(iy)+0.5)/double(m_ny)*2.0-1.0);
}

OrbitCalculator::~OrbitCalculator()
{
  writeSlice(m_nslice);
}

bool OrbitCalculator::writeSlice(unsigned new_islice)
{
  if(m_islice == new_islice)return false;
  if(m_islice != m_nslice)
    {
      std::ostringstream fn; fn << "sl" << m_islice << ".dat";
      std::cerr << "Writing: " << fn.str() << '\n';
      std::ofstream str(fn.str().c_str());
      for(unsigned iy=0;iy<m_ny;iy++)
	{
	  for(unsigned ix=0;ix<m_nx;ix++)
	    {
	      if(ix)str << ' ';
	      str << m_data[iel(ix,iy)];
	    }
	  str << '\n';
	}

      m_data.clear();
    }
  if(new_islice<m_nslice)
    {
      m_data.resize(m_nx*m_ny,0);
      m_islice = new_islice;
      return true;
    }
  else
    {
      m_islice = m_nslice;
      return false;
    }
}

void OrbitCalculator::visitElement(unsigned irow, FT2& t)
{
  if(m_t0<=0.0)m_t0=t.t_start;
  if(t.t_start<m_t0)return;
  unsigned islice = unsigned(floor((t.t_start-m_t0)/m_tslice));
  bool newslice = writeSlice(islice);
  if(islice>=m_nslice)return;
  m_islice=islice;
  
  if(t.in_saa)return;

  if(newslice)
    {
      m_j2o = Vec3D(0,0,-t.pole_ra);
      m_j2o &= Vec3D(0,-(M_PI_2-t.pole_dec),0);

      Vec3D galcent = Vec3D::makeLatLon(0, 0);
      galcent.rotate(m_g2j);
      galcent.rotate(m_j2o);

      double theta = atan2(galcent.y(),galcent.x());
      m_j2o &= Vec3D(0,0,-theta);

      galcent = Vec3D::makeLatLon(0, 0);
      galcent.rotate(m_g2j);
      galcent.rotate(m_j2o);

      Vec3D galpole = Vec3D::makeLatLon(M_PI_2, 0);
      galpole.rotate(m_g2j);
      galpole.rotate(m_j2o);
  
      std::cout 
	<< m_islice << ' ' 
	<< r2d(galcent.lat()) << ' ' << r2d(galcent.lon()) << ' '
	<< r2d(galpole.lat()) << ' ' << r2d(galpole.lon()) << ' '
	<< m_j2o.x() << ' ' << m_j2o.y() << ' ' << m_j2o.z() << '\n';
    }

  Vec3D scz  = Vec3D::makeLatLon(t.scz_dec, t.scz_ra);
  scz.rotate(m_j2o);

  for(unsigned ix=0;ix<m_nx;ix++)
    for(unsigned iy=0;iy<m_ny;iy++)
      {
	Vec3D xy = Vec3D::makeLatLon(m_y[iy],m_x[ix]);
	double costheta = xy*scz;
	//Vec3D cp = xy^scz;
	//double sintheta = cp.norm();
	//theta = 
	if(costheta>=0.2)m_data[iel(ix,iy)] += ea(costheta);
      }
}

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] ft2_file"
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);

  std::string command_line(*argv);
  for(int iarg=1;iarg<argc;iarg++)
    {
      command_line+=" ";
      command_line+=argv[iarg];
    }

  // --------------------------------------------------------------------------
  // PROCESS OPTIONS
  // --------------------------------------------------------------------------

  VSOptions options(argc, argv, true);
  
  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  // --------------------------------------------------------------------------
  // FINISH OPTIONS PROCESSING
  // --------------------------------------------------------------------------

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
	if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  int arg_req = 3;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string ft2_fn(*argv);
  argc--, argv++;

  double tslice(30); VSDataConverter::fromString(tslice,*argv);
  argc--, argv++;
  if(tslice<1.0)tslice=1.0;

  unsigned nslice(1); VSDataConverter::fromString(nslice,*argv);
  argc--, argv++;

  double t0(0);
  if(argc)
    {
      VSDataConverter::fromString(t0,*argv);
      argc--, argv++;
    }

  // --------------------------------------------------------------------------
  // PROCESS THE EVENTS
  // --------------------------------------------------------------------------

  FITSVectorVisitor<FT2>* visitor(0);
  visitor = new OrbitCalculator(t0,tslice,nslice);

  //  FITSVectorDispatcher<FT2> dispatcher(visitor);
  //  unsigned nraw = dispatcher.dispatchVector(ft2_fn, FT2::tableName());

  Magic7Dispatcher dispatcher(visitor);
  dispatcher.dispatchVector(ft2_fn);
}
