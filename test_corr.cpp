//-*-mode:c++; mode:font-lock;-*-

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cassert>

#include "VERITAS/RandomNumbers.hpp"
#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"

#include "MyMinuit.hpp"
#include "Accumulator.hpp"

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

inline double SQR(const double& x) { return x*x; }

class CorrelationWithDataError: public Optimizable
{
public:
  CorrelationWithDataError(const std::vector<double>& x,
			   const std::vector<double>& dx,
			   const std::vector<double>& y,
			   const std::vector<double>& dy,
			   bool force_zero=false);
  virtual ~CorrelationWithDataError();
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);

private:
  bool                m_force_zero;
  unsigned            m_N;
  std::vector<double> m_x;
  std::vector<double> m_dx;
  std::vector<double> m_y;
  std::vector<double> m_dy;
};

CorrelationWithDataError::
CorrelationWithDataError(const std::vector<double>& x,
			 const std::vector<double>& dx,
			 const std::vector<double>& y,
			 const std::vector<double>& dy,
			 bool force_zero):
  Optimizable(), m_force_zero(force_zero), 
  m_N(x.size()), m_x(x), m_dx(dx), m_y(y), m_dy(dy)
{
  assert(m_dx.size() == m_N);
  assert(m_N  == m_N);
  assert(m_dy.size() == m_N);
}

CorrelationWithDataError::~CorrelationWithDataError()
{
  // nothing to see here
}

bool CorrelationWithDataError::canCalcDFDP()
{
  return true;//false;
}

double CorrelationWithDataError::f(const std::vector<double>& p)
{
  const double mx(p[0]);
  const double my(p[1]);
  const double sx(p[2]);
  const double sy(p[3]);
  const double r(m_force_zero ? 0 : p[4]);
  Accumulator logL;
  for(unsigned i=0; i<m_N; i++)
    {
      const double xi(m_x[i]);
      const double dxi(m_dx[i]);
      const double yi(m_y[i]);
      const double dyi(m_dy[i]);

      const double Sxxi(SQR(sx)+SQR(dxi));
      const double Syyi(SQR(sy)+SQR(dyi));
      const double Sxyi(r*sx*sy);

      const double detSi = Sxxi*Syyi-Sxyi*Sxyi;

      const double Mxxi = Syyi/detSi;
      const double Myyi = Sxxi/detSi;
      const double Mxyi = -Sxyi/detSi;

      const double vxi = xi-mx;
      const double vyi = yi-my;

      logL.add(-0.5*std::log(detSi));
      logL.add(-0.5*(Mxxi*SQR(vxi) + Myyi*SQR(vyi) + 2*Mxyi*vxi*vyi));
    }
  return -logL.sum();
}

double CorrelationWithDataError::
dfdp(const std::vector<double>& p, unsigned iparam)
{
  const double mx(p[0]);
  const double my(p[1]);
  const double sx(p[2]);
  const double sy(p[3]);
  const double r(m_force_zero ? 0 : p[4]);
  Accumulator dLogL_dX;
  for(unsigned i=0; i<m_N; i++)
    {
      const double xi(m_x[i]);
      const double dxi(m_dx[i]);
      const double yi(m_y[i]);
      const double dyi(m_dy[i]);

      const double Sxxi(SQR(sx)+SQR(dxi));
      const double Syyi(SQR(sy)+SQR(dyi));
      const double Sxyi(r*sx*sy);

      const double detSi = Sxxi*Syyi-Sxyi*Sxyi;

      const double Mxxi = Syyi/detSi;
      const double Myyi = Sxxi/detSi;
      const double Mxyi = -Sxyi/detSi;

      const double vxi = xi-mx;
      const double vyi = yi-my;

      switch(iparam)
	{
	case 0:
	  dLogL_dX.add(Mxxi*vxi + Mxyi*vyi);
	  break;
	case 1:
	  dLogL_dX.add(Myyi*vyi + Mxyi*vxi);
	  break;
	case 2:
	  {
	    const double dSxxi_dsx(2.0*sx);
	    const double dSxyi_dsx(r*sy);

	    const double ddetSi_dsx = dSxxi_dsx*Syyi-2.0*Sxyi*dSxyi_dsx;

	    const double dMxxi_dsx = -Syyi/SQR(detSi)*ddetSi_dsx;
	    const double dMyyi_dsx = dSxxi_dsx/detSi-Sxxi/SQR(detSi)*ddetSi_dsx;
	    const double dMxyi_dsx = -dSxyi_dsx/detSi+Sxyi/SQR(detSi)*ddetSi_dsx;
	    dLogL_dX.add(-0.5*ddetSi_dsx/detSi);
	    dLogL_dX.add(-0.5*(dMxxi_dsx*SQR(vxi) + dMyyi_dsx*SQR(vyi) 
			       + 2*dMxyi_dsx*vxi*vyi));
	  }
	  break;
	case 3:
	  {
	    const double dSyyi_dsy(2.0*sy);
	    const double dSxyi_dsy(r*sx);

	    const double ddetSi_dsy = Sxxi*dSyyi_dsy-2.0*Sxyi*dSxyi_dsy;

	    const double dMxxi_dsy = dSyyi_dsy/detSi-Syyi/SQR(detSi)*ddetSi_dsy;
	    const double dMyyi_dsy = -Sxxi/SQR(detSi)*ddetSi_dsy;
	    const double dMxyi_dsy = -dSxyi_dsy/detSi+Sxyi/SQR(detSi)*ddetSi_dsy;

	    dLogL_dX.add(-0.5*ddetSi_dsy/detSi);
	    dLogL_dX.add(-0.5*(dMxxi_dsy*SQR(vxi) + dMyyi_dsy*SQR(vyi)
			       + 2*dMxyi_dsy*vxi*vyi));
	  }
	  break;
	case 4:
	  {
	    const double dSxyi_dr(sx*sy);

	    const double ddetSi_dr = -2.0*Sxyi*dSxyi_dr;

	    const double dMxxi_dr = -Syyi/SQR(detSi)*ddetSi_dr;
	    const double dMyyi_dr = -Sxxi/SQR(detSi)*ddetSi_dr;
	    const double dMxyi_dr = -dSxyi_dr/detSi + Sxyi/SQR(detSi)*ddetSi_dr;
	    
	    dLogL_dX.add(-0.5*ddetSi_dr/detSi);
	    dLogL_dX.add(-0.5*(dMxxi_dr*SQR(vxi) + dMyyi_dr*SQR(vyi) 
			       + 2*dMxyi_dr*vxi*vyi));
	  }
	  break;
	}
    }
  return -dLogL_dX.sum();
}

unsigned CorrelationWithDataError::numParam()
{
  if(m_force_zero)return 4;
  else return 5;
}

std::string CorrelationWithDataError::pName(unsigned iparam)
{
  switch(iparam)
    {
    case 0:
      return "x_mean";
    case 1:
      return "y_mean";
    case 2:
      return "x_sigma";
    case 3:
      return "y_sigma";
    case 4:
      return "rho";
    default:
      return std::string();
    }
}

double CorrelationWithDataError::p0(unsigned iparam)
{
  switch(iparam)
    {
    case 0:
      {
	double x(0);
	for(unsigned i=0;i<m_N;i++)x += m_x[i];
	return x/double(m_N);
      }
    case 1:
      {
	double y(0);
	for(unsigned i=0;i<m_N;i++)y += m_y[i];
	return y/double(m_N);
      }      
    case 2:
      {
	double x(0);
	double xx(0);
	for(unsigned i=0;i<m_N;i++)
	  x += m_x[i], xx += SQR(m_x[i]);
	return std::sqrt(xx/double(m_N) - SQR(x/double(m_N)));
      }
    case 3:
      {
	double y(0);
	double yy(0);
	for(unsigned i=0;i<m_N;i++)
	  y += m_y[i], yy += SQR(m_y[i]);
	return std::sqrt(yy/double(m_N) - SQR(y/double(m_N)));
      }
    default:
      return 0;
    }  
}

double CorrelationWithDataError::plo(unsigned iparam)
{
  switch(iparam)
    {
    case 0:
      {
	double xmin(m_x[0]);
	for(unsigned i=1;i<m_N;i++)xmin=std::min(xmin,m_x[i]);
	return xmin;
      }
    case 1:
      {
	double ymin(m_y[0]);
	for(unsigned i=1;i<m_N;i++)ymin=std::min(ymin,m_y[i]);
	return ymin;
      }
    case 2:
    case 3:
      return 0;
    case 4:
      return -1;
    default:
      return 0;
    }
}

double CorrelationWithDataError::phi(unsigned iparam)
{
  switch(iparam)
    {
    case 0:
      {
	double xmax(m_x[0]);
	for(unsigned i=1;i<m_N;i++)xmax=std::max(xmax,m_x[i]);
	return xmax;
      }
    case 1:
      {
	double ymax(m_y[0]);
	for(unsigned i=1;i<m_N;i++)ymax=std::max(ymax,m_y[i]);
	return ymax;
      }      
    case 2:
      {
	double xmin(m_x[0]);
	double xmax(m_x[0]);
	for(unsigned i=1;i<m_N;i++)
	  xmin=std::min(xmin,m_x[i]), xmax=std::max(xmax,m_x[i]);
	return xmax-xmin;
      }
    case 3:
      {
	double ymin(m_y[0]);
	double ymax(m_y[0]);
	for(unsigned i=1;i<m_N;i++)
	  ymin=std::min(ymin,m_y[i]), ymax=std::max(ymax,m_y[i]);
	return ymax-ymin;
      }
    case 4:
      return 1.0;
    default:
      return 0;
    }  
}

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] data_file"
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

  int verbose = 0;
  if(options.find("v", "Verbose.")!=VSOptions::FS_NOT_FOUND)
    verbose = 1;
  if(options.find("vv", "Very verbose.")!=VSOptions::FS_NOT_FOUND)
    verbose = 2;

  bool interactive = false;
  if(options.find("interactive", "Run minimizer in interactive mode.")
     !=VSOptions::FS_NOT_FOUND)
    interactive=true;

  unsigned nrun = 0;
  if(options.findWithValue("nrun", nrun,
			   "Set the number or runs to make either shuffling "
			   "(default) the y-data randomly or doing a full "
			   "simulated dataset (if sim option is given.")
     ==VSOptions::FS_FOUND_BUT_WITHOUT_VALUE)
    nrun=1;

  bool sim = false;
  std::vector<double> defsimparam(9);
  defsimparam[0] = 0;
  defsimparam[1] = 0;
  defsimparam[2] = 1;
  defsimparam[3] = 1;
  defsimparam[4] = 0;
  defsimparam[5] = 1;
  defsimparam[6] = 1;
  defsimparam[7] = 0;
  defsimparam[8] = 0;
  std::vector<double> simparam(defsimparam);
  if(options.findWithValue("sim", simparam,
			   "Enable simulation mode and set the values of the "
			   "7 simulation parameters. These are xmean, "
			   "ymean, xsig, ysig, xycorr, xmeaserr, ymeaserr, "
			   "xmeanmeaserr, ymeanmeaserr.")
     !=VSOptions::FS_NOT_FOUND)
    sim = true;

  while(sim && simparam.size()<defsimparam.size())
    simparam.push_back(defsimparam[simparam.size()]);

  unsigned nwrite = 0;
  if(options.findWithValue("writesim", nwrite,
			   "In simulation mode, do not base number of data "
			   "points on entries in file, rather generate given "
			   "numbers and write them to file.")
     ==VSOptions::FS_FOUND_BUT_WITHOUT_VALUE)
    nwrite=1;
  
  if(options.find("interactive", "Run minimizer in interactive mode.")
     !=VSOptions::FS_NOT_FOUND)
    interactive=true;

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

  int arg_req = 1;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string data_fn(*argv);
  argc--, argv++;

  std::vector<double> x;
  std::vector<double> dx;
  std::vector<double> y;
  std::vector<double> dy;

  if(sim && nwrite>0)
    {
      x.resize(nwrite);
      y.resize(nwrite);
      dx.resize(nwrite);
      dy.resize(nwrite);
    }
  else
    {
      std::ifstream str(data_fn.c_str());
      std::string line;
      std::getline(str,line);
      while(str)
	{
	  std::istringstream lstr(line);
	  double _x;
	  double _dx;
	  double _y;
	  double _dy;
	  lstr >> _x >> _dx >> _y >> _dy;
	  x.push_back(_x);
	  dx.push_back(_dx);
	  y.push_back(_y);
	  dy.push_back(_dy);
	  std::getline(str,line);
	}
    }

  RandomNumbers* rng(0);
  if(sim || nrun)
    rng = new RandomNumbers(RandomNumbers::defaultFilename());

  unsigned irun=0;
  do
    {
      if(sim)
	{
	  unsigned nxy(y.size());
	  for(unsigned ixy=0;ixy<nxy;ixy++)
	    {
	      double mx = rng->Normal()*simparam[2];
	      double my = rng->Normal()*simparam[3] + mx*simparam[4];
	      mx += simparam[0];
	      my += simparam[1];
	      double ex = 0.5*(1.0 + std::abs(rng->Normal()))*simparam[5] +
		std::sqrt(std::abs(mx))*simparam[7];
	      double ey = 0.5*(1.0 + std::abs(rng->Normal()))*simparam[6] +
		std::sqrt(std::abs(my))*simparam[8];
	      x[ixy] = mx + ex*rng->Normal();
	      y[ixy] = my + ey*rng->Normal();
	      dx[ixy] = ex;
	      dy[ixy] = ey;
	    }

	  if(nwrite>0 && irun==0)
	    {
	      std::ofstream str(data_fn.c_str());
	      for(unsigned ixy=0;ixy<nxy;ixy++)
		str << x[ixy] << ' ' << dx[ixy] << ' '
		    << y[ixy] << ' ' << dy[ixy] << '\n';
	    }
	}
      else if(nrun)
	{
	  unsigned ny(y.size());
	  for(unsigned iy=0;iy<(ny-1);iy++)
	    {
	      unsigned jy = (rng->UInt32()%(ny-iy)); // Small bias here
	      std::swap(y[iy],y[iy+jy]);
	      std::swap(dy[iy],dy[iy+jy]);
	    }
	}
      
      double logL_ML;
      double logL_ZH;

      if(1)
	{
	  CorrelationWithDataError cde(x,dx,y,dy,false);
	  MyMinuit M(&cde,verbose-1);

	  try
	    {
	      M.minimize(interactive);
	      logL_ML = M.fVal();
	    }
	  catch(const std::string& s)
	    {
	      std::cerr << s << '\n';
	      logL_ML = 0;
	    }

	  if(nrun<=1)
	    std::cout
	      << "ML hypothesis:\n"
	      << " X_mean: " << M.pVal(0) << " +/- " << M.pErr(0) << '\n'
	      << " Y_mean: " << M.pVal(1) << " +/- " << M.pErr(1) << '\n'
	      << " X_sigma: " << M.pVal(2) << " +/- " << M.pErr(2) << '\n'
	      << " Y_sigma: " << M.pVal(3) << " +/- " << M.pErr(3) << '\n'
	      << " Corr: " << M.pVal(4) << " +/- " << M.pErr(4) << '\n'
	      << " Likelihood: " << logL_ML << "\n\n";
	  else
	    std::cout
	      << M.pVal(0) << ' ' << M.pErr(0) << ' '
	      << M.pVal(1) << ' ' << M.pErr(1) << ' '
	      << M.pVal(2) << ' ' << M.pErr(2) << ' '
	      << M.pVal(3) << ' ' << M.pErr(3) << ' '
	      << M.pVal(4) << ' ' << M.pErr(4) << ' ';
	}
      
      if(1)
	{
	  CorrelationWithDataError cde(x,dx,y,dy,true);
	  MyMinuit M(&cde,verbose-1);

	  try
	    {
	      M.minimize(interactive);
	      logL_ZH = M.fVal();
	    }
	  catch(const std::string& s)
	    {
	      std::cerr << s << '\n';
	      logL_ZH = 0;
	    }
	  
	  if(nrun<=1)
	    std::cout
	      << "NULL hypothesis:\n"
	      << " X_mean: " << M.pVal(0) << " +/- " << M.pErr(0) << '\n'
	      << " Y_mean: " << M.pVal(1) << " +/- " << M.pErr(1) << '\n'
	      << " X_sigma: " << M.pVal(2) << " +/- " << M.pErr(2) << '\n'
	      << " Y_sigma: " << M.pVal(3) << " +/- " << M.pErr(3) << '\n'
	      << " Likelihood: " << logL_ZH << "\n\n";
	  else
	    std::cout
	      << M.pVal(0) << ' ' << M.pErr(0) << ' '
	      << M.pVal(1) << ' ' << M.pErr(1) << ' '
	      << M.pVal(2) << ' ' << M.pErr(2) << ' '
	      << M.pVal(3) << ' ' << M.pErr(3) << ' ';
	}
      
      irun++;

      if(nrun<=1)std::cout << "2DeltaLogL: ";
      else std::cout << logL_ML << ' ' << logL_ZH << ' ';
      std::cout << 2.0*(logL_ZH-logL_ML) << '\n';
    }while(irun<nrun);

  delete rng;
}
