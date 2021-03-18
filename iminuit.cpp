//-*-mode:c++; mode:font-lock;-*-

#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>

#include "VERITAS/VSOptions.hpp"

#include "MyMinuit.hpp"

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname
         << " [options] command_argument"
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

class PopenOptimizable: public Optimizable
{
public:
  typedef quad<std::string,double,double,double> PDef;

  PopenOptimizable(int argc, char** argv,
		   const std::vector<PDef>& pdef, int verbose);
  virtual ~PopenOptimizable();
  virtual unsigned numParam();
  virtual std::string pName(unsigned iparam);
  virtual bool canCalcDFDP();
  virtual double f(const std::vector<double>& p);
  virtual double dfdp(const std::vector<double>& p, unsigned iparam);
  virtual double p0(unsigned iparam);
  virtual double plo(unsigned iparam);
  virtual double phi(unsigned iparam);
private:
  std::string lookupVariable(const std::string& variable,
			     const std::vector<double>& p);

  std::string m_arg;
  std::vector<PDef> m_pdef;
  int m_verbose;
};

PopenOptimizable::
PopenOptimizable(int argc, char** argv, const std::vector<PDef>& pdef,
		 int verbose):
  Optimizable(), m_arg(), m_pdef(pdef), m_verbose(verbose)
{
  if(argc==0)
    throw std::string("What you talkin' 'bout Willis?");
  std::ostringstream stream;
  stream << *argv;
  argv++,argc--;
  while(argc)
    {   
      stream << ' ' << *argv;
      argv++,argc--;
    }
  m_arg = stream.str();
}

PopenOptimizable::~PopenOptimizable()
{
  // nothing to see here
}

bool PopenOptimizable::canCalcDFDP()
{
  return false;
}

double PopenOptimizable::f(const std::vector<double>& p)
{
  std::string command;

  bool escape = false;
  bool dollar_first = false;
  bool dollar = false;
  bool bracket = false;
  std::string variable;
  for(std::string::const_iterator ichar = m_arg.begin();
      ichar != m_arg.end(); ichar++)
    {
      if(escape)
	{
	  command += *ichar;
	  escape = false;
	  continue;
	}

      if(dollar_first)
	{
	  dollar_first = false;
	  if(*ichar=='(')
	    {
	      bracket=true;
	      continue;
	    }
	  else 
	    dollar=true;
	}

      if(bracket)
	{
	  if(*ichar==')')
	    {
	      bracket=false;
	      command += lookupVariable(variable,p);
	      variable.clear();
	      continue;
	    }
	  else
	    {
	      variable += *ichar;
	      continue;
	    }
	}

      if(dollar)
	{
	  if(isalnum(*ichar)||(*ichar=='_'))
	    {
	      variable += *ichar;
	      continue;
	    }
	  else
	    {
	      command += lookupVariable(variable,p);
	      variable.clear();
	      dollar=false;
	    }
	}

      if(*ichar == '\\')
	{
	  escape = true;
	}
      else if(*ichar == '$')
	{
	  dollar_first = true;
	}
      else
	{
	  command += *ichar;
	}
    }

  if(dollar)
    {
      command += lookupVariable(variable,p);
      variable.clear();
    }

  if(m_verbose>1)std::cout << command << std::flush;
  FILE* fp = popen(command.c_str(), "r");
  if(fp==0)throw std::string("Could not open subprocess");
  double x = 0;
  assert(fscanf(fp,"%lf",&x));
  pclose(fp);
  if(m_verbose>1)
    std::cout << " = " << VSDataConverter::toString(x) << std::endl
	      << std::flush;;
  
  return x;
}

double PopenOptimizable::dfdp(const std::vector<double>& p, unsigned iparam)
{
  throw std::string("What you talkin' 'bout Willis?");
}

unsigned PopenOptimizable::numParam()
{
  return m_pdef.size();
}

std::string PopenOptimizable::pName(unsigned iparam)
{
  return m_pdef[iparam].first;
}

double PopenOptimizable::p0(unsigned iparam)
{
  return m_pdef[iparam].second;
}

double PopenOptimizable::plo(unsigned iparam)
{
  return m_pdef[iparam].third;
}

double PopenOptimizable::phi(unsigned iparam)
{
  return m_pdef[iparam].fourth;
}

std::string PopenOptimizable::lookupVariable(const std::string& variable,
					     const std::vector<double>& p)
{
  unsigned np = m_pdef.size();
  for(unsigned ip=0;ip<np;ip++)
    {
      if(variable == m_pdef[ip].first)
	return VSDataConverter::toString(p[ip]);
    }
  return "";
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

  int verbose = 1;
  if(options.find("q","Print no messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 0;
  if(options.find("v","Print really verbose messages during analysis.")
     !=VSOptions::FS_NOT_FOUND)
    verbose = 2;
  
  std::vector<PopenOptimizable::PDef> p;
  options.findWithValue("p",p,"Define parameters. Should be given as a comma "
			"separated list of four entries: "
			"\"pname,p0,plo,phi\", where \"pname\" defines the "
			"name of the parameter, p0 the initial guess and plo "
			"& phi the low and high  bound for the search range.");

  bool interactive = false;
  options.findBoolValue("interactive", interactive, true, 
		       "Run MINUIT interactively.");
  options.findBoolValue("i", interactive, true, 
		       "Run MINUIT interactively.");

  double tol=1;
  options.findWithValue("tol",tol,
			"Set the MINUIT tolerance in parts per thousand.");

  double err=0.5;
  options.findWithValue("err",err,
			"Set the MINUIT error definition.");

  unsigned strategy=1;
  options.findWithValue("strategy",strategy,
			"Set the MINUIT minimization strategy.");
  options.findWithValue("str",strategy,
			"Set the MINUIT minimization strategy.");

  double epsilon=0.0;
  options.findWithValue("epsilon",epsilon,
			"Set the MINUIT floating point precision.");
  options.findWithValue("eps",epsilon,
			"Set the MINUIT floating point precision.");  

  std::list<std::string> minos;
  if(options.findWithValue("minos",minos,
			   "Calculate the asymmetric errors for some "
			   "parameters using MINOS. Parameters should be given "
			   "as a comma separated list of names. If no "
			   "parameter is listed then MINOS will be run on all "
			   "parameters.")
     == VSOptions::FS_FOUND_BUT_WITHOUT_VALUE) minos.push_back("*");

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

  if(argc < 1)
    {
      std::cerr << progname << ": need at least one argument." << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }
  
  // --------------------------------------------------------------------------
  // CRANK UP MINUIT
  // --------------------------------------------------------------------------

  PopenOptimizable fn(argc,argv,p,verbose);

  if(verbose)
    {
      std::cout << "Parameters:\n";
      unsigned np = fn.numParam();
      for(unsigned ip=0; ip<np; ip++)
	std::cout << "- " << ip+1 << ' '
		  << fn.pName(ip) << ' ' << fn.p0(ip) << ' '
		  << fn.plo(ip) << ' ' << fn.phi(ip) << '\n';
    }
  
  try
    {
      MyMinuit minuit(&fn,verbose==0?-1:verbose);
      if(epsilon>0)minuit.setEpsilon(epsilon);
      minuit.minimize(interactive,tol,err,"",false,strategy);

      while(!minos.empty())
	{
	  std::string param(minos.front());
	  minos.pop_front();
	  if(param == "*")
	    {
	      for(unsigned ip=p.size();ip;ip--)
		minos.push_front(p[ip-1].first);
	    }
	  else
	    {
	      unsigned ip = 0;
	      while(ip<p.size())
		{
		  if(param == p[ip].first)break;
		  ip++;
		}
	      if(ip == p.size())
		{
		  std::cerr << "Minos: unknown parameter: " << param << '\n';
		}
	      else
		{
		  double perr_minus;
		  double perr_plus;
		  minuit.computeMinosErrors(ip, perr_minus, perr_plus);
		  std::cout << "Minos: " << param << ' ' 
			    << perr_minus << ' ' << perr_plus << '\n';
		}
	    }
	}
    }
  catch(const std::string& x)
    {
      std::cerr << "Caught exception: " << x << '\n';
    }
}
