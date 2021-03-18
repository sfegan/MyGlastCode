//-*-mode:c++; mode:font-lock;-*-

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "tip/TipException.h"

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"

#include "FT1.hpp"
#include "FT2.hpp"
#include "Util.hpp"
#include "Analysis.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

class SOBinWriter: public SOEventVisitor
{
public:
  SOBinWriter(const std::string& fn);
  ~SOBinWriter();
  virtual void visitSOEvent(unsigned irow, FT1& event, Vec3D& so);
private:
  int m_fd;
};

#define W(x,y,z) assert(write(x,y,z)==z)

SOBinWriter::SOBinWriter(const std::string& fn): SOEventVisitor(), m_fd()
{
  m_fd = open(fn.c_str(), O_WRONLY|O_CREAT|O_TRUNC, 
	      S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH);
  if(m_fd<0)throw std::string("Could not open "+fn);
  unsigned magic(0xFEEDBEEF);
  W(m_fd,&magic,sizeof(magic));
}

SOBinWriter::~SOBinWriter()
{
  close(m_fd);
}

void SOBinWriter::visitSOEvent(unsigned irow, FT1& event, Vec3D& so)
{
  W(m_fd,&event.run_id,sizeof(event.run_id));
  W(m_fd,&event.event_id,sizeof(event.event_id));
  double d;
  d = so.x(); W(m_fd,&d,sizeof(d));
  d = so.y(); W(m_fd,&d,sizeof(d));
  d = so.z(); W(m_fd,&d,sizeof(d));
}

class SOPrinter: public SOEventVisitor
{
public:
  SOPrinter(std::ostream& str): SOEventVisitor(), m_str(str) { }
  virtual ~SOPrinter();
  virtual void visitSOEvent(unsigned irow, FT1& event, Vec3D& so);
private:
  std::ostream& m_str;
};
  
SOPrinter::~SOPrinter()
{
  // nothing to see here
}

void SOPrinter::visitSOEvent(unsigned irow, FT1& event, Vec3D& so)
{
#if 1
  m_str
    << event.run_id << ' ' << event.event_id << ' '
    << VSDataConverter::toString(so.x()) << ' '
    << VSDataConverter::toString(so.y()) << ' '
    << VSDataConverter::toString(so.z()) << '\n';
#endif

#if 0
  Vec3D d = Vec3D::makePolar(event.theta,event.phi);
  d.rotate(so);
  std::cout << r2d(event.ra) << ' ' << r2d(event.dec) << ' '
	    << r2d(d.phi()) << ' ' << r2d(M_PI_2-d.theta()) << ' '
	    << r2d(sphere_dist(event.ra, event.dec, d.phi(), M_PI_2-d.theta()))
	    << '\n';
#endif

#if 0
  Vec3D d = Vec3D::makePolar(M_PI_2-event.dec,event.ra);
  d.rotate(-so);
  
  std::cout << r2d(event.ra) << ' ' << r2d(event.dec) << ' '
	    << r2d(d.phi()) << ' ' << r2d(M_PI_2-d.theta()) << ' '
	    << r2d(sphere_dist(event.ra, event.dec, d.phi(), M_PI_2-d.theta()))
	    << '\n';
#endif
}

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] ft1_file ft2_file"
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);

  // --------------------------------------------------------------------------
  // PROCESS OPTIONS
  // --------------------------------------------------------------------------

  VSOptions options(argc, argv, true);
  
  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool binary = false;
  if(options.find("binary","Use binary output file.")!=VSOptions::FS_NOT_FOUND)
    binary=true;

  std::string ofn("");
  options.findWithValue("o",ofn,"Set output file name.");

  if(binary && ofn.empty())
    {
      std::cerr << "Must supply output file name in binary mode.\n";
      exit(EXIT_FAILURE);
    }

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

  int arg_req = 2;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string ft1_fn(*argv);
  argc--, argv++;


  std::string ft2_fn(*argv);
  argc--, argv++;

  // --------------------------------------------------------------------------
  // LOAD THE SPACECRAFT (FT2) FILE
  // --------------------------------------------------------------------------

  std::vector<FT2> ft2_vec;
  FT2::loadFromFITS(ft2_vec, ft2_fn);
#if 0 
  unsigned nft2 = 
  std::cout << "Loaded " << nft2 << " spacecraft entries from "
	    << ft2_fn << '\n';
#endif

  // --------------------------------------------------------------------------
  // PROCESS THE EVENTS
  // --------------------------------------------------------------------------

  FITSVectorVisitor<FT1>* visitor(0);

  SOEventVisitor* printer(0);
  std::ofstream* streamf(0);

  if(binary)
    printer = new SOBinWriter(ofn);
  else if(ofn.empty() || ofn=="STDOUT" || ofn=="-")
    printer = new SOPrinter(std::cout);
  else
    {
      streamf = new std::ofstream(ofn.c_str());
      printer = new SOPrinter(*streamf);
    }
  
  QuadInterpSOCalculator socalc(printer, ft2_vec);
  visitor = &socalc;

  FITSVectorDispatcher<FT1> dispatcher(visitor);
  unsigned nrawevent = dispatcher.dispatchVector(ft1_fn, FT1::tableName());
  unsigned nsoevent = socalc.nSOEventsDispatched();

  std::cerr << "Oriented " << nsoevent << '/' << nrawevent << " events\n";

  delete printer;
  delete streamf;
}
