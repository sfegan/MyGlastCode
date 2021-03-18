#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <map>
#include <list>
#include <vector>
#include <cstdlib>

//#define IMINUIT_USE_READLINE
#ifdef IMINUIT_USE_READLINE
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include "MyMinuit.hpp"

ParameterizableFunction::~ParameterizableFunction()
{
  // nothing to see here
}

std::string ParameterizableFunction::pName(unsigned iparam)
{
  std::ostringstream pline;
  pline << "P" << iparam+1;
  return pline.str();
}

ProbDensity::~ProbDensity()
{
  // nothing to see here
}

IntProbDensity::~IntProbDensity()
{
  // nothing to see here
}

Optimizable::~Optimizable()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
//
// PARAMETER OPTIMIZER
//
// ----------------------------------------------------------------------------

ParameterModifier::
ParameterModifier(Optimizable* fn): 
  Optimizable(), m_fn(fn), m_pmap(), m_pfixed(), m_p0(), m_phi(), m_plo()
{
  unsigned nparam = fn->numParam();
  for(unsigned iparam=0;iparam<nparam; iparam++)
    m_pmap[fn->pName(iparam)]=iparam;
}

ParameterModifier::~ParameterModifier()
{
  // nothing to see here
}

unsigned ParameterModifier::numParam()
{
  return m_pmap.size()-m_pfixed.size();
}

unsigned ParameterModifier::iFnParam(unsigned iparam) const
{
  unsigned ifnparam = iparam;
  for(std::map<unsigned, double>::const_iterator ifixed = m_pfixed.begin();
      ifixed != m_pfixed.end(); ifixed++)
    if(ifixed->first<=ifnparam)ifnparam++;
  return ifnparam;
}

std::string ParameterModifier::pName(unsigned iparam)
{
  return m_fn->pName(iFnParam(iparam));
}

bool ParameterModifier::canCalcDFDP()
{
  return m_fn->canCalcDFDP();
}

double ParameterModifier::f(const std::vector<double>& p)
{
  std::vector<double> fnp(m_pmap.size());
  for(unsigned iparam=0;iparam<p.size();iparam++)
    fnp[iFnParam(iparam)] = p[iparam];
  for(std::map<unsigned, double>::const_iterator p_fix = m_pfixed.begin();
      p_fix != m_pfixed.end(); p_fix++)
    fnp[p_fix->first] = p_fix->second;
  return m_fn->f(fnp);
}

double ParameterModifier::dfdp(const std::vector<double>& p, unsigned iparam)
{
  std::vector<double> fnp(m_pmap.size());
  for(unsigned jparam=0;jparam<p.size();jparam++)
    fnp[iFnParam(jparam)] = p[jparam];
  for(std::map<unsigned, double>::const_iterator p_fix = m_pfixed.begin();
      p_fix != m_pfixed.end(); p_fix++)
    fnp[p_fix->first] = p_fix->second;
  return m_fn->dfdp(fnp,iFnParam(iparam));
}

double ParameterModifier::p0(unsigned iparam)
{
  unsigned ifnparam = iFnParam(iparam);
  std::map<unsigned, double>::const_iterator p_find = m_p0.find(ifnparam);
  if(p_find == m_p0.end())return p_find->second;
  else return m_fn->p0(iparam);
}

double ParameterModifier::plo(unsigned iparam)
{
  unsigned ifnparam = iFnParam(iparam);
  std::map<unsigned, double>::const_iterator p_find = m_plo.find(ifnparam);
  if(p_find == m_plo.end())return p_find->second;
  else return m_fn->plo(iparam);
}

double ParameterModifier::phi(unsigned iparam)
{
  unsigned ifnparam = iFnParam(iparam);
  std::map<unsigned, double>::const_iterator p_find = m_phi.find(ifnparam);
  if(p_find == m_phi.end())return p_find->second;
  else return m_fn->phi(iparam);
}

// ----------------------------------------------------------------------------
//
// MY MINUIT
//
// ----------------------------------------------------------------------------

MyMinuit::MyMinuit(Optimizable* fn, int verbose)
  : m_fn(fn), m_verbose(verbose)
{ 
  const int i5=5, i6=6, i7=7;
  mninit_(&i5, &i6, &i7);
}


std::vector<std::string> MyMinuit::pNames() const
{
  unsigned np(m_fn->numParam());
  std::vector<std::string> pv(np);
  for(unsigned ip=0;ip<np;ip++)pv[ip] = m_fn->pName(ip);
  return pv;
}

std::string MyMinuit::pName(unsigned iparam) const
{
  return m_fn->pName(iparam);
}

void MyMinuit::setEpsilon(double epsilon)
{
  std::ostringstream streaml;
  streaml << "SET EPS " << epsilon;
  doCmd(streaml.str());
}

void MyMinuit::
minimize(bool interactive, double tol, double err,
	 const std::string& command_str, bool nohesse, unsigned strategy)
{
  int error_flag;
  
  if(m_verbose > 0)
    {
      const int i5=5, i6=6, i7=7;
      mintio_(&i5, &i6, &i7);
    }

  std::ostringstream pline;
  if(!interactive || m_verbose>=0)
    {
      pline << "SET PRINT " << m_verbose;
      doCmd(pline.str()); // Set verbosity of Minuit
    }
  doCmd("SET NOWARN");
  pline.str(std::string());
  pline << "SET ERR " << err;
  doCmd(pline.str());
  if(m_fn->canCalcDFDP())
    doCmd("SET GRAD 1");  // Use gradient calculated by fcn

  unsigned nparam = m_fn->numParam();
  for(unsigned iparam=0; iparam<nparam; iparam++)
    {
      std::string name = m_fn->pName(iparam);
      int iparam_fortran = iparam+1;
      double p0 = m_fn->p0(iparam);
      double plo = m_fn->plo(iparam);
      double phi = m_fn->phi(iparam);
      double scale = 1.0;
      mnparm_(&iparam_fortran, name.c_str(), &p0, &scale, &plo, &phi,
	      &error_flag, name.size());
    }
  
  tol *= 1000.0;
  if(interactive)minimization_program_interactive(tol,command_str);
  else minimization_program_fixed(tol,nohesse,strategy);
}

void MyMinuit::
minimization_program_fixed(double tol, bool nohesse, unsigned strategy)
{
  int error_flag;
  std::ostringstream streaml;
  streaml << "SET STR " << strategy;
  error_flag = doCmd(streaml.str());
  streaml.str(std::string());
  streaml << "MIN 0 " << tol;
  error_flag = doCmd(streaml.str());
  if (error_flag == 4) 
    throw std::string("Minuit abnormal termination. (No convergence?)");
  else if (error_flag > 0)
    throw std::string("Minuit bad command line");
  if(!nohesse)doCmd("HESSE 0");
}

void MyMinuit::
minimization_program_interactive(double tol, const std::string& command_str)
{
  std::vector<std::string> commands;
  if(!command_str.empty())
    {
      // Parse command text into individual command strings
      std::string word;
      bool escape = false;
      for(std::string::const_iterator ichar = command_str.begin();
	  ichar != command_str.end(); ichar++)
	{
	  if(escape)
	    {
	      if(*ichar != ',')word += '\\'; // save escape char for later
	      word += *ichar;
	      escape = false;
	      //std::cout << *ichar << " : " << word << '\n';
	      continue;
	    }

	  if(*ichar == '\\')
	    {
	      escape = true;
	    }
	  else if(*ichar == ',')
	    {
	      if(!word.empty())commands.push_back(word);
	      word.clear();
	    }
	  else
	    {
	      word += *ichar;
	    }
	  //std::cout << *ichar << " : " << word << '\n';
	}

      if(!word.empty())commands.push_back(word);
    }

  std::vector<std::string>::const_iterator icommand(commands.begin());

  unsigned nparam = numParam();
  std::map<std::string,unsigned> pnames;
  for(unsigned iparam=0;iparam<nparam;iparam++)
    pnames[pName(iparam)] = iparam+1;

  if(m_verbose)
    std::cout 
      << std::endl
      << "iMINUIT: Interactive MINUIT optimizer. Type \".h\" for help."
      << std::endl
      << std::endl;

  // Read in commands and send them to Minuit
  std::istream* input(&std::cin);
  std::list<std::istream*> input_list;
  bool go = true;
  int last_cmd_status = 0;

  while(go)
    {
      std::string command;
      if(icommand != commands.end())
	{
	  // First process all commands given through constructor
	  command = *icommand;
	  icommand++;
	}
      else
	{
	  // Then process commands from stdin or open files
	  if(input->eof())
	    {
	      if(input != &std::cin)
		{
		  delete input;
		  input = 0;
		}

	      if(input_list.empty())
		{
		  go = false;
		  continue;
		}
	      else
		{
		  input = input_list.front();
		  input_list.pop_back();
		}
	    }

	  if(input == &std::cin)
	    {
#ifdef IMINUIT_USE_READLINE
	      char* line = readline("iMINUIT> ");
	      if(!line)
		{
		  if(input_list.empty())
		    {
		      go = false;
		      continue;
		    }
		  else
		    {
		      input = input_list.front();
		      input_list.pop_back();
		    }
		}
	      if(*line)add_history(line);
	      command = line;
	      free(line);
#else
	      std::cout << "iMINUIT> " << std::flush;
	      std::getline(*input,command);
#endif
	    }
	  else
	    {
	      std::cout << "iMINUIT> " << std::flush;
	      std::getline(*input,command);
	      std::cout << command << std::endl;
	    }
	}

      std::vector<std::string> cmdwords = 
	splitCommand(command, pnames, tol);

    reruncommand:
      if(cmdwords.empty())continue;
      unsigned ncmdword = cmdwords.size();

      if(cmdwords[0][0]=='.')
	{
	  switch(cmdwords[0].size())
	    {
	    case 1:
	      continue;
	    case 2:
	      break;
	    default:
	      std::cout 
		<< "iMINUIT: unrecognized command: " << cmdwords[0]
		<< std::endl;
	      continue;
	    }

	  switch(cmdwords[0][1])
	    {
	    case 'h':
	      printHelp(true);
	      break;
	    case '?':
	      printHelp(false);
	      break;
	    case 'l':
	      printParameters();
	      break;
	    case 'e':
	      for(unsigned iword=1;iword<ncmdword;iword++)
		{
		  if(iword!=1)std::cout << ' ';
		  std::cout << cmdwords[iword];
		}
	      std::cout << std::endl;
	      break;
	    case 'q':
	    case 'x':
	      go = false;
	      break;
	    case 'X':
	      if(ncmdword!=1)
		throw std::string("MyMinuit exception requested");
	      else
		{
		  std::ostringstream cstream;
		  cstream << cmdwords[1];
		  for(unsigned iword=2; iword<ncmdword; iword++)
		    cstream << ' ' << cmdwords[iword];
		  throw cstream.str();
		}
	      break;
	    case 'i':
	      if(ncmdword != 2)
		{
		  std::cout
		    << "iMINUIT: error: .i command needs FILENAME"
		    << std::endl;
		  continue;
		}
	      input_list.push_back(input);
	      input = new std::ifstream(cmdwords[1].c_str());
	      if(!input->good())
		{
		  std::cout
		    << "iMINUIT: could not open file: " << cmdwords[1]
		    << std::endl;
		  delete input;
		  input = input_list.front();
		  input_list.pop_back();
		}
	      break;
	    case 'r':
	      if(ncmdword != 2)
		{
		  std::cout
		    << "iMINUIT: error: .r command needs FILENAME"
		    << std::endl;
		  continue;
		}
	      input_list.push_back(input);
	      input = new std::ifstream(cmdwords[1].c_str());
	      if(!input->good())
		{
		  std::cout
		    << "iMINUIT: could not open file: " << cmdwords[1]
		    << std::endl;
		  delete input;
		  input = input_list.front();
		  input_list.pop_back();
		}
	      else
		{
		  std::istream* str = input_list.front();
		  input_list.pop_back();
		  if(str != &std::cin)delete str;
		}
	      break;
	    case 'D':
	    case 'd':
	      if(ncmdword==1)dumpState(std::cout, cmdwords[0][1]=='D');
	      else
		{
		  std::ofstream of(cmdwords[1].c_str()); 		    
		  if(of.good())
		    {
		      dumpState(of,cmdwords[0][1]=='D');
		      std::cout << "iMINUIT: wrote state to file: "
				<< cmdwords[1] << std::endl;
		    }
		  else
		    std::cout << "iMINUIT: could not open file: "
			      << cmdwords[1] << std::endl;
		}
	      break;
	    case 'c':
	      if(ncmdword < 3)
		{
		  std::cout
		    << "iMINUIT: error: .c command needs X and Y variables"
		    << std::endl;
		  continue;
		}
	      else
		{
		  std::istringstream s;
		  std::vector<double> x;
		  std::vector<double> y;
		  unsigned ix = 1;
		  unsigned iy = 2;
		  unsigned N = 20;
		  std::istringstream(cmdwords[1]) >> ix;
		  std::istringstream(cmdwords[2]) >> iy;
		  ix--;
		  iy--;
		  if(ncmdword == 4)std::istringstream(cmdwords[3]) >> N;
		  if(computeMinosContour(ix, iy, x, y, N)>0)
		    for(unsigned ip=0;ip<x.size();ip++)
		      std::cout << "iMINUIT_MC: " 
				<< x[ip] << ' ' << y[ip] << '\n';
		}
	      break;
	    case 's':
	      if(last_cmd_status != 4)
		{
		  cmdwords.erase(cmdwords.begin());
		  goto reruncommand;
		}
	    case 'S':
	      if(last_cmd_status == 4)
		{
		  cmdwords.erase(cmdwords.begin());
		  goto reruncommand;
		}
	    default:
	      std::cout 
		<< "iMINUIT: unrecognized command: " << cmdwords[0]
		<< std::endl;
	    }
	}
      else
	{
	  std::ostringstream cstream;
	  cstream << cmdwords[0];
	  for(unsigned iword=1; iword<ncmdword; iword++)
	    cstream << ' ' << cmdwords[iword];
	  last_cmd_status = doCmd(cstream.str());
	}
    }

  while(!input_list.empty())
    {
      delete input;
      input = input_list.front();
      input_list.pop_front();
    }
  if(input != &std::cin)
    {
      delete input;
      input = 0;
    }
}

double MyMinuit::fVal() const
{
  MinStat ms(getMinimizationStatus());
  return ms.fmin;
}

std::vector<double> MyMinuit::pVals() const
{
  std::vector<ValAndErr> ve_vec(getValuesAndErrors());
  unsigned nparam(ve_vec.size());
  std::vector<double> p_vec(nparam);
  for(unsigned iparam=0;iparam<nparam;iparam++)
    p_vec[iparam] = ve_vec[iparam].val;
  return p_vec;
}

double MyMinuit::pVal(unsigned iparam) const
{
  ValAndErr ve(getValueAndError(iparam));
  return ve.val;
}

std::vector<double> MyMinuit::pErrs() const
{
  std::vector<ValAndErr> ve_vec(getValuesAndErrors());
  unsigned nparam(ve_vec.size());
  std::vector<double> perr_vec(nparam);
  for(unsigned iparam=0;iparam<nparam;iparam++)
    perr_vec[iparam] = ve_vec[iparam].err;
  return perr_vec;
}

double MyMinuit::pErr(unsigned iparam) const
{
  ValAndErr ve(getValueAndError(iparam));
  return ve.err;
}

std::vector<std::vector<double> > MyMinuit::pCovMtx() const
{
  int nparam = m_fn->numParam();
  std::vector<double> elements(nparam*nparam);
  mnemat_(&elements.front(), &nparam);
  std::vector<std::vector<double> >
    matrix(nparam,std::vector<double>(nparam));
  for(int irow(0), iel(0); irow < nparam; irow++)
    for (int icol(0); icol < nparam; icol++)
      matrix[irow][icol] = elements[iel++];
  return matrix;
}

MyMinuit::MinStat MyMinuit::getMinimizationStatus() const
{
  MinStat ms;
  mnstat_(&ms.fmin, &ms.fedm, &ms.errdef, &ms.npari, &ms.nparx, &ms.istat);
  return ms;
}

std::vector<MyMinuit::ValAndErr>
MyMinuit::getValuesAndErrors() const
{
  MinStat ms = getMinimizationStatus();
  unsigned nparam = ms.nparx;
  std::vector<MyMinuit::ValAndErr> ve_vec(nparam);
  for(unsigned iparam=0;iparam<nparam;iparam++)
    {
      char buffer[11];
      buffer[10]='\0';
      MyMinuit::ValAndErr& ve(ve_vec[iparam]);
      int j = iparam+1;
      int intvar;
      mnpout_(&j, buffer, &ve.val, &ve.err, &ve.bnd_lo, &ve.bnd_hi, &intvar,
	      sizeof(buffer)/sizeof(*buffer)-1);
      ve.pname = pName(iparam);
      ve.mname = buffer;
      mnerrs_(&j, &ve.err_plus, &ve.err_minus, &ve.err_parabolic,
		&ve.global_corr);
    }
  return ve_vec;
}

MyMinuit::ValAndErr MyMinuit::getValueAndError(unsigned iparam) const
{
  char buffer[11];
  buffer[10]='\0';
  MyMinuit::ValAndErr ve;
  int j = iparam+1;
  int intvar;
  mnpout_(&j, buffer, &ve.val, &ve.err, &ve.bnd_lo, &ve.bnd_hi, &intvar,
	  sizeof(buffer)/sizeof(*buffer)-1);
  ve.pname = pName(iparam);
  ve.mname = buffer;
  mnerrs_(&j, &ve.err_plus, &ve.err_minus, &ve.err_parabolic,
	  &ve.global_corr);
  return ve;
}

int MyMinuit::
computeMinosErrors(unsigned iparam, double& perr_minus, double& perr_plus)
{
  std::ostringstream s;
  s << "MINOS 0 " << iparam+1;
  int err_flag = doCmd(s.str());
  MyMinuit::ValAndErr ve = getValueAndError(iparam);
  perr_minus = ve.err_minus;
  perr_plus  = ve.err_plus;
  return err_flag;
}

int MyMinuit::
computeMinosContour(unsigned xparam, unsigned yparam,
		    std::vector<double>& x, std::vector<double>& y,
		    unsigned n)
{
  int par1 = xparam+1;
  int par2 = yparam+1;
  x.resize(n);
  y.resize(n);
  int npt_req = n;
  int npt_fnd = 0;
  mncont_(&s_fcn, &par1, &par2, &npt_req, 
	  &x.front(), &y.front(), &npt_fnd,
	  static_cast<void*>(this));
  if(npt_fnd>0)
    {
      x.resize(npt_fnd);
      y.resize(npt_fnd);
    }
  else
    { 
      x.clear();
      y.clear();
    }
  return npt_fnd;
}

// ----------------------------------------------------------------------------
//
// PRIVATE FUNCTIONS
//
// ----------------------------------------------------------------------------

int MyMinuit::doCmd(const std::string& command)
{
  // Pass a command string to Minuit
  int errorFlag = 0;  
  mncomd_(&s_fcn, command.c_str(), &errorFlag, static_cast<void*>(this),
	  command.length());
  return errorFlag;
}

void MyMinuit::fcn(int* npar, double* grad, double* fcnval,
		   double* pval, int* iflag)
{
  std::vector<double> p(pval, pval+m_fn->numParam());
  *fcnval = m_fn->f(p);
  if (*iflag == 2)for(int i=0; i<*npar; i++)grad[i] = m_fn->dfdp(p,i);
}

void MyMinuit::
s_fcn(int* npar, double* grad, double* fcnval,
      double* pval, int* iflag, void* futil)
{
  // This is the function that Minuit minimizes.
  // Minuit thinks futil is a function pointer. It's been hijacked to
  // be a pointer to MyMinuit, so this static function can use it.
  MyMinuit* mm = static_cast<MyMinuit*>(futil);
  mm->fcn(npar,grad,fcnval,pval,iflag);
}

std::vector<std::string> MyMinuit::
splitCommand(const std::string& cmdtxt,
	     const std::map<std::string,unsigned>& parameters,
	     double tol) const
{
  std::vector<std::string> words;
  std::string word;
  bool escape = false;
  for(std::string::const_iterator ichar = cmdtxt.begin();
      ichar != cmdtxt.end(); ichar++)
    {
      if(escape)
	{
	  word += *ichar;
	  escape = false;
	  continue;
	}

      if(*ichar == '#')
	{
	  break;
	}
      else if(*ichar == '\\')
	{
	  escape = true;
	}
      else if(isspace(*ichar))
	{
	  if(!word.empty())
	    {
	      words.push_back(word);
	      word.clear();
	    }
	}
      else
	{
	  word += *ichar;
	}
    }
  if(!word.empty())words.push_back(word);
  for(unsigned iword=0;iword<words.size();iword++)
    if(words[iword][0]=='$')
      {
	std::string var = words[iword].substr(1);
	std::ostringstream vstream;
	if(var == "TOL" || var == "tol")
	  vstream << tol;
	else if(var.substr(0,5) == "chi2(" || var.substr(0,5) == "CHI2(")
	  {
	    std::istringstream stream(var.substr(5));
	    double P = 0.6827;
	    double n = 1.0;
	    char c;
	    stream >> P >> c >> n;
	    vstream << 0.5*invchisquaredistribution(n,1.0-P);
	  }
	else
	  {
	    std::map<std::string,unsigned>::const_iterator iparam =
	      parameters.find(var);
	    if(iparam != parameters.end())
	      vstream << iparam->second;
	    else
	      {
		std::cout << "iMINUIT: variable not found $" << var
			  << std::endl;
		continue;
	      }
	  }
	words[iword]=vstream.str();
      }
  return words;
}

void MyMinuit::printHelp(bool long_help) const
{
  if(long_help)
    std::cout << "iMINUIT: Interactive Minuit optimizer\n\
\n\
This optimizer allows the user to run CERN's MINUIT package\n\
interactively or to specify a set of commands used by the minimizer\n\
in advance (i.e. to customize the MINUIT \"program\" with respect to\n\
what is used by the standard ScienceTools Minuit optimizer).\n\
\n\
Commands starting with a period are interpreted here, while all\n\
others are given directly to MINUIT. See the CERN manual for MINUIT\n\
for more details of what commands it accepts, or type \"HELP\" to\n\
access MINUIT's interactive help system.\n\
\n\
The following commands are recognized here:\n\
\n\
.h     - print this help message.\n\
\n\
.?     - print a short command summary.\n";
    else std::cout << "\
.h     - print full help message.\n\
\n\
.?     - print this command summary.\n";
      std::cout << "\
\n\
.i XXX - read command from file (name given as XXX) and execute them as\n\
         if they were typed here. Continue interactively after command in\n\
         file are completed.\n\
\n\
.r XXX - read command from file (name given as XXX) and execute them as\n\
         if they were typed here. Terminate after complete.\n\
\n\
.e XXX - echo everything following the \".e\". Useful for expanding \"$\"\n\
         variables (see explanation below).\n\
\n\
.l     - list parameter names and numbers being used by MINUIT.\n\
\n\
.d XXX - dump state of minimization to file XXX (or to standard output if\n\
.D XXX   no name is given). Use \".d\" if you do not want the covariance\n\
         matrix, or \".D\" if you do.\n\
\n\
.c X Y - compute a 2-D error contour using MINOS. X and Y are the variables\n\
         to use on the X and Y axes of the contour respectively.\n\
\n\
.x/.q  - terminate minimization session and allow ScienceTools to\n\
         complete its work (such as calculating TS values in gtlike).\n\
         This should only be done when you are happy you have\n\
         optimized the model and calculated the errors.\n\
\n\
.X XXX - throw an exception with message \"XXX\".\n\
\n\
.s XXX - execute the command \"XXX\" only if the last Minuit command ended\n\
.S XXX   in success (\".s\") or in failure (\".S\").\n";
  if(long_help)
    std::cout << "\
\n\
In addition to these commands iMINUIT will replace some words starting\n\
with a \"$\". \"$TOL\" is translated into the tolerance requested by the\n\
caller. \"$CHI2(P,N)\" is translated into the difference in a Gaussian\n\
likelihood function that corresponds to a probability level of P with N\n\
degrees of freedom. Secondly, \"$XXX\", where XXX is the name of one of the\n\
parameters (listed by \".l\") is replaced with the number of the parameter\n\
that MINUIT will use.\n\
\n\
For example, to emulate command that the noninteractive MINUIT optimizer\n\
issues, type the following commands:\n\
\n\
iMINUIT> MIN 0 0.001\n\
iMINUIT> HESSE\n\
iMINUIT> .q\n\
\n\
A more complex example, using MINOS to estimate the errors of two\n\
parameters, and plotting a contour of the 68% region defined by the two\n\
of them is given below:\n\
\n\
iMINUIT> SET STRATEGY 2       # use strictest convergence and error calc\n\
iMINUIT> MIN 0 0.001OL        # run MIGRAD with no call limit\n\
iMINUIT> HESSE                # calculate error matrix\n\
iMINUIT> .l                   # get list of parameters\n\
iMINUIT> MINOS 0 4 5          # run MINOS on parameters 4&5, with no call limit\n\
iMINUIT> .D fit.dat           # dump parameter values, errors and corr mat\n\
iMINUIT> SET ERR $CHI2(.68,2) # set error level to 68% with 2 DOF\n\
iMINUIT> SET PRINT 1          # set verbosity so contour points are printed\n\
iMINUIT> MNC 4 5              # calculate confidence contour on parms. 4&5\n\
iMINUIT> .q                   # return to ScienceTools\n\
\
" << std::endl;
}

void MyMinuit::printParameters() const
{
  std::vector<std::string> pnames = pNames();
  std::cout << "List of parameters:" << std::endl;
  for (unsigned iparam = 0; iparam<pnames.size(); iparam++)
    std::cout << iparam+1 << ' ' << pnames[iparam] << std::endl;
}

//! Print dump of current state
void MyMinuit::dumpState(std::ostream& str, 
			 bool dumpCovarianceMatrix) const
{
  MinStat ms = getMinimizationStatus();

  {
    std::ostringstream lstr;
    int np = 3;
    if(ms.fmin>0)np = unsigned(ceil(std::log10(ms.fmin)))+3;
    if(np<3)np=3;
    lstr << std::setprecision(np) << ms.fmin;
    str << "iMINUIT_MS " 
	<< lstr.str() << ' '
	<< ms.fedm << ' ' << ms.errdef << ' '
	<< ms.npari << ' ' << ms.nparx << ' ' << ms.istat << '\n';
  }
    
  std::vector<ValAndErr> ve = getValuesAndErrors();
  unsigned nve = ve.size();
    
  unsigned maxusize = 0;
  unsigned maxpsize = 0;
  for(unsigned ive=0; ive<nve; ive++)
    if(ve[ive].pname.size()>maxpsize)maxpsize=ve[ive].pname.size();

  unsigned np = 4;
  {
    std::ostringstream lstr;
    lstr << "iMINUIT_CC " << "# " << ' '
	 << std::setw(maxusize) << std::left << "P_NAME" << ' '
      // << std::setw(maxpsize) << std::left << "P_NAME" << ' '
      // << std::setw(10)"  << std::left << "P_NAME" << ' '
	 << std::setw(np+7) << std::left << " LIM_LO" << ' '
	 << std::setw(np+7) << std::left << " LIM_HI" << ' '
	 << std::setw(np+7) << std::left << " VAL" << ' '
	 << std::setw(np+7) << std::left << " ERR" << ' '
	 << std::setw(np+7) << std::left << " ERR_PAR" << ' '
	 << std::setw(np+7) << std::left << " ERR_MINUS" << ' '
	 << std::setw(np+7) << std::left << " ERR_PLUS" << ' '
	 << std::setw(6) << std::left << "CORR" 
	 << '\n';
    str << lstr.str();
  }
    
  for(unsigned ive=0; ive<nve; ive++)
    {
      std::ostringstream lstr;
      lstr << "iMINUIT_VE " << std::setw(2) << std::left << ive+1 << ' '
	   << std::setw(maxpsize) << std::left << ve[ive].pname << ' ' 
	// << ve[ive].mname << ' '
	   << std::scientific << std::showpos
	   << std::setprecision(np) << ve[ive].bnd_lo << ' ' 
	   << std::setprecision(np) << ve[ive].bnd_hi << ' '
	   << std::setprecision(np) << ve[ive].val << ' ' 
	   << std::setprecision(np) << ve[ive].err << ' ' 
	   << std::setprecision(np) << ve[ive].err_parabolic << ' '
	   << std::setprecision(np) << ve[ive].err_minus << ' '
	   << std::setprecision(np) << ve[ive].err_plus << ' ' 
	   << std::fixed << std::noshowpos << ve[ive].global_corr
	   << '\n';
	str << lstr.str();
    }

  if(dumpCovarianceMatrix)
    {
      np = 6;
      std::vector<std::vector<double> > matrix = 
	MyMinuit::pCovMtx();
      for(unsigned irow=0; irow<nve; irow++)
	{
	  str << "iMINUIT_CM";
	  for(unsigned icol=0; icol<nve; icol++)
	    {
	      std::ostringstream lstr;
	      lstr 
		<< ' '
		<< std::scientific << std::showpos << std::setprecision(np)
		<< matrix[irow][icol];
	      str << lstr.str();
	    }
	  str << '\n';
	}
    }

  str << std::flush;
}

/*************************************************************************
**************************************************************************
**************************************************************************
***                                                                    ***
*** ALL REMAINING CODE IN THIS FILE ARE TO CALCULATE CHI2 PROBABILITY  ***
***                                                                    ***
*** http://www.alglib.net/specialfunctions/distributions/chisquare.php ***
***                                                                    ***
*** Code has following copyright and licencing information             ***
***                                                                    ***
**************************************************************************
**************************************************************************
*************************************************************************/


/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
      pseudocode.

See subroutines comments for additional copyrights.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

	   
/*************************************************************************
Inverse of complemented Chi-square distribution

Finds the Chi-square argument x such that the integral
from x to infinity of the Chi-square density is equal
to the given cumulative probability y.

This is accomplished using the inverse gamma integral
function and the relation

   x/2 = igami( df/2, y );

ACCURACY:

See inverse incomplete gamma function


Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/

double MyMinuit::invchisquaredistribution(double v, double y)
{
  double result;
  if(y<0 || y>1 || v<1)
    throw std::string("Domain error in InvChiSquareDistribution");
  result = 2*invincompletegammac(0.5*v, y);
  return result;
}

/*************************************************************************
Inverse of complemented imcomplete gamma integral

Given p, the function finds x such that

 igamc( a, x ) = p.

Starting with the approximate value

        3
 x = a t

 where

 t = 1 - d - ndtri(p) sqrt(d)

and

 d = 1/9a,

the routine performs up to 10 Newton iterations to find the
root of igamc(a,x) - p = 0.

ACCURACY:

Tested at random a, p in the intervals indicated.

               a        p                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
   IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
   IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/

double MyMinuit::invincompletegammac(double a, double y0)
{
  double result;
  double igammaepsilon;
  double iinvgammabignumber;
  double x0;
  double x1;
  double x;
  double yl;
  double yh;
  double y;
  double d;
  double lgm;
  double dithresh;
  int i;
  int dir;
  double tmp;

  igammaepsilon = 0.000000000000001;
  iinvgammabignumber = 4503599627370496.0;
  x0 = iinvgammabignumber;
  yl = 0;
  x1 = 0;
  yh = 1;
  dithresh = 5*igammaepsilon;
  d = 1/(9*a);
  y = 1-d-invnormaldistribution(y0)*sqrt(d);
  x = a*y*y*y;
  lgm = lngamma(a, tmp);
  i = 0;
  while(i<10)
    {
      if( x>x0||x<x1 )
	{
	  d = 0.0625;
	  break;
	}
      y = incompletegammac(a, x);
      if( y<yl||y>yh )
	{
	  d = 0.0625;
	  break;
	}
      if( y<y0 )
	{
	  x0 = x;
	  yl = y;
	}
      else
	{
	  x1 = x;
	  yh = y;
	  }
      d = (a-1)*std::log(x)-x-lgm;
      if( d<-709.78271289338399 )
	{
	  d = 0.0625;
	  break;
	}
      d = -std::exp(d);
      d = (y-y0)/d;
      if( fabs(d/x)<igammaepsilon )
	{
	  result = x;
	  return result;
	}
      x = x-d;
      i = i+1;
    }
  if( x0==iinvgammabignumber )
    {
      if( x<=0 )
	{
	  x = 1;
	}
      while(x0==iinvgammabignumber)
	{
	  x = (1+d)*x;
	  y = incompletegammac(a, x);
	  if( y<y0 )
	    {
	      x0 = x;
	      yl = y;
	      break;
	    }
	  d = d+d;
	}
    }
  d = 0.5;
  dir = 0;
  i = 0;
  while(i<400)
    {
      x = x1+d*(x0-x1);
      y = incompletegammac(a, x);
      lgm = (x0-x1)/(x1+x0);
      if( fabs(lgm)<dithresh )
	{
	  break;
	}
      lgm = (y-y0)/y0;
      if( fabs(lgm)<dithresh )
	{
	  break;
	}
      if( x<=0.0 )
	{
	  break;
	}
      if( y>=y0 )
	{
	  x1 = x;
	  yh = y;
	  if( dir<0 )
	    {
	      dir = 0;
	      d = 0.5;
	    }
	  else
	    {
	      if( dir>1 )
		{
		  d = 0.5*d+0.5;
		}
	      else
		{
		  d = (y0-yl)/(yh-yl);
		}
	    }
	  dir = dir+1;
	}
      else
	{
	  x0 = x;
	  yl = y;
	  if( dir>0 )
	    {
	      dir = 0;
	      d = 0.5;
	    }
	  else
	    {
	      if( dir<-1 )
		{
		  d = 0.5*d;
		}
	      else
		{
		  d = (y0-yl)/(yh-yl);
		}
	    }
	  dir = dir-1;
	}
      i = i+1;
    }
  result = x;
  return result;
}

/*************************************************************************
Complemented incomplete gamma integral

The function is defined by


 igamc(a,x)   =   1 - igam(a,x)

                           inf.
                             -
                    1       | |  -t  a-1
              =   -----     |   e   t   dt.
                   -      | |
                  | (a)    -
                            x


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

Tested at random a, x.
               a         x                      Relative error:
arithmetic   domain   domain     # trials      peak         rms
   IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
   IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/

double MyMinuit::incompletegammac(double a, double x)
{
  double result;
  double igammaepsilon;
  double igammabignumber;
  double igammabignumberinv;
  double ans;
  double ax;
  double c;
  double yc;
  double r;
  double t;
  double y;
  double z;
  double pk;
  double pkm1;
  double pkm2;
  double qk;
  double qkm1;
  double qkm2;
  double tmp;

  igammaepsilon = 0.000000000000001;
  igammabignumber = 4503599627370496.0;
  igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
  if( x<=0||a<=0 )
    {
      result = 1;
      return result;
    }
  if( x<1||x<a )
    {
      result = 1-incompletegamma(a, x);
      return result;
    }
  ax = a*std::log(x)-x-lngamma(a, tmp);
  if( ax<-709.78271289338399 )
    {
      result = 0;
      return result;
    }
  ax = std::exp(ax);
  y = 1-a;
  z = x+y+1;
  c = 0;
  pkm2 = 1;
  qkm2 = x;
  pkm1 = x+1;
  qkm1 = z*x;
  ans = pkm1/qkm1;
  do
    {
      c = c+1;
      y = y+1;
      z = z+2;
      yc = y*c;
      pk = pkm1*z-pkm2*yc;
      qk = qkm1*z-qkm2*yc;
      if( qk!=0 )
	{
	  r = pk/qk;
	  t = fabs((ans-r)/r);
	  ans = r;
	}
      else
	{
	  t = 1;
	}
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if( fabs(pk)>igammabignumber )
	{
	  pkm2 = pkm2*igammabignumberinv;
	  pkm1 = pkm1*igammabignumberinv;
	  qkm2 = qkm2*igammabignumberinv;
	  qkm1 = qkm1*igammabignumberinv;
	}
    }
  while(t>igammaepsilon);
  result = ans*ax;
  return result;
}

/*************************************************************************
Incomplete gamma integral

The function is defined by

                          x
                           -
                  1       | |  -t  a-1
 igam(a,x)  =   -----     |   e   t   dt.
                 -      | |
                | (a)    -
                          0


In this implementation both arguments must be positive.
The integral is evaluated by either a power series or
continued fraction expansion, depending on the relative
values of a and x.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30       200000       3.6e-14     2.9e-15
   IEEE      0,100      300000       9.9e-14     1.5e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/

double MyMinuit::incompletegamma(double a, double x)
{
  double result;
  double igammaepsilon;
  double ans;
  double ax;
  double c;
  double r;
  double tmp;
    
  igammaepsilon = 0.000000000000001;
  if( x<=0||a<=0 )
    {
      result = 0;
      return result;
    }
  if( x>1&&x>a )
    {
      result = 1-incompletegammac(a, x);
      return result;
    }
  ax = a*std::log(x)-x-lngamma(a, tmp);
  if( ax<-709.78271289338399 )
    {
      result = 0;
      return result;
    }
  ax = std::exp(ax);
  r = a;
  c = 1;
  ans = 1;
  do
    {
      r = r+1;
      c = c*x/r;
      ans = ans+c;
    }
  while(c/ans>igammaepsilon);
  result = ans*ax/a;
  return result;
}

/*************************************************************************
Natural logarithm of gamma function

Input parameters:
    X       -   argument

Result:
    logarithm of the absolute value of the Gamma(X).

Output parameters:
    SgnGam  -   sign(Gamma(X))

Domain:
    0 < X < 2.55e305
    -2.55e305 < X < 0, X is not an integer.

ACCURACY:
arithmetic      domain        # trials     peak         rms
   IEEE    0, 3                 28000     5.4e-16     1.1e-16
   IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
The error criterion was relative when the function magnitude
was greater than one but absolute when it was less than one.

The following test used the relative error criterion, though
at certain points the relative error could be much higher than
indicated.
   IEEE    -200, -4             10000     4.8e-16     1.3e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
*************************************************************************/

double MyMinuit::lngamma(double x, double& sgngam)
{
  double result;
  double a;
  double b;
  double c;
  double p;
  double q;
  double u;
  double w;
  double z;
  int i;
  double logpi;
  double ls2pi;
  double tmp;

  sgngam = 1;
  logpi = 1.14472988584940017414;
  ls2pi = 0.91893853320467274178;
  if( x<-34.0 )
    {
      q = -x;
      w = lngamma(q, tmp);
      p = int(floor(q));
      i = int(floor(p+0.5));
      if( i%2==0 )
	{
	  sgngam = -1;
	}
      else
	{
	  sgngam = 1;
	}
      z = q-p;
      if( z>0.5 )
	{
	  p = p+1;
	  z = p-q;
	}
      z = q*std::sin(M_PI*z);
      result = logpi-std::log(z)-w;
      return result;
    }
  if( x<13 )
    {
      z = 1;
      p = 0;
      u = x;
      while(u>=3)
	{
	  p = p-1;
	  u = x+p;
	  z = z*u;
	}
      while(u<2)
	{
	  z = z/u;
	  p = p+1;
	  u = x+p;
	}
      if( z<0 )
	{
	  sgngam = -1;
	  z = -z;
	}
      else
	{
	  sgngam = 1;
	}
      if( u==2 )
	{
	  result = std::log(z);
	  return result;
	}
      p = p-2;
      x = x+p;
      b = -1378.25152569120859100;
      b = -38801.6315134637840924+x*b;
      b = -331612.992738871184744+x*b;
      b = -1162370.97492762307383+x*b;
      b = -1721737.00820839662146+x*b;
      b = -853555.664245765465627+x*b;
      c = 1;
      c = -351.815701436523470549+x*c;
      c = -17064.2106651881159223+x*c;
      c = -220528.590553854454839+x*c;
      c = -1139334.44367982507207+x*c;
      c = -2532523.07177582951285+x*c;
      c = -2018891.41433532773231+x*c;
      p = x*b/c;
      result = std::log(z)+p;
      return result;
    }
  q = (x-0.5)*std::log(x)-x+ls2pi;
  if( x>100000000 )
    {
      result = q;
      return result;
    }
  p = 1/(x*x);
  if( x>=1000.0 )
    {
      q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
    }
  else
    {
      a = 8.11614167470508450300*0.0001;
      a = -5.95061904284301438324*0.0001+p*a;
      a = 7.93650340457716943945*0.0001+p*a;
      a = -2.77777777730099687205*0.001+p*a;
      a = 8.33333333333331927722*0.01+p*a;
      q = q+a/x;
    }
  result = q;
  return result;
}

/*************************************************************************
Inverse of Normal distribution function

Returns the argument, x, for which the area under the
Gaussian probability density function (integrated from
minus infinity to x) is equal to y.


For small arguments 0 < y < exp(-2), the program computes
z = sqrt( -2.0 * log(y) );  then the approximation is
x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
There are two rational functions P/Q, one for 0 < y < exp(-32)
and the other for y up to exp(-2).  For larger arguments,
w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).

ACCURACY:

                     Relative error:
arithmetic   domain        # trials      peak         rms
   IEEE     0.125, 1        20000       7.2e-16     1.3e-16
   IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*************************************************************************/

double MyMinuit::invnormaldistribution(double y0)
{
  double result;
  double expm2;
  double s2pi;
  double x;
  double y;
  double z;
  double y2;
  double x0;
  double x1;
  int code;
  double p0;
  double q0;
  double p1;
  double q1;
  double p2;
  double q2;
  
  expm2 = 0.13533528323661269189;
  s2pi = 2.50662827463100050242;
  if( y0<=0 )
    {
      result = -DBL_MAX;
      return result;
    }
  if( y0>=1 )
    {
      result = DBL_MAX;
      return result;
    }
  code = 1;
  y = y0;
  if( y>1.0-expm2 )
    {
      y = 1.0-y;
      code = 0;
    }
  if( y>expm2 )
    {
      y = y-0.5;
      y2 = y*y;
      p0 = -59.9633501014107895267;
      p0 = 98.0010754185999661536+y2*p0;
      p0 = -56.6762857469070293439+y2*p0;
      p0 = 13.9312609387279679503+y2*p0;
      p0 = -1.23916583867381258016+y2*p0;
      q0 = 1;
      q0 = 1.95448858338141759834+y2*q0;
      q0 = 4.67627912898881538453+y2*q0;
      q0 = 86.3602421390890590575+y2*q0;
      q0 = -225.462687854119370527+y2*q0;
      q0 = 200.260212380060660359+y2*q0;
      q0 = -82.0372256168333339912+y2*q0;
      q0 = 15.9056225126211695515+y2*q0;
      q0 = -1.18331621121330003142+y2*q0;
      x = y+y*y2*p0/q0;
      x = x*s2pi;
      result = x;
      return result;
    }
  x = std::sqrt(-2.0*std::log(y));
  x0 = x-std::log(x)/x;
  z = 1.0/x;
  if( x<8.0 )
    {
      p1 = 4.05544892305962419923;
      p1 = 31.5251094599893866154+z*p1;
      p1 = 57.1628192246421288162+z*p1;
      p1 = 44.0805073893200834700+z*p1;
      p1 = 14.6849561928858024014+z*p1;
      p1 = 2.18663306850790267539+z*p1;
      p1 = -1.40256079171354495875*0.1+z*p1;
      p1 = -3.50424626827848203418*0.01+z*p1;
      p1 = -8.57456785154685413611*0.0001+z*p1;
      q1 = 1;
      q1 = 15.7799883256466749731+z*q1;
      q1 = 45.3907635128879210584+z*q1;
      q1 = 41.3172038254672030440+z*q1;
      q1 = 15.0425385692907503408+z*q1;
      q1 = 2.50464946208309415979+z*q1;
      q1 = -1.42182922854787788574*0.1+z*q1;
      q1 = -3.80806407691578277194*0.01+z*q1;
      q1 = -9.33259480895457427372*0.0001+z*q1;
      x1 = z*p1/q1;
    }
  else
    {
      p2 = 3.23774891776946035970;
      p2 = 6.91522889068984211695+z*p2;
      p2 = 3.93881025292474443415+z*p2;
      p2 = 1.33303460815807542389+z*p2;
      p2 = 2.01485389549179081538*0.1+z*p2;
      p2 = 1.23716634817820021358*0.01+z*p2;
      p2 = 3.01581553508235416007*0.0001+z*p2;
      p2 = 2.65806974686737550832*0.000001+z*p2;
      p2 = 6.23974539184983293730*0.000000001+z*p2;
      q2 = 1;
      q2 = 6.02427039364742014255+z*q2;
      q2 = 3.67983563856160859403+z*q2;
      q2 = 1.37702099489081330271+z*q2;
      q2 = 2.16236993594496635890*0.1+z*q2;
      q2 = 1.34204006088543189037*0.01+z*q2;
      q2 = 3.28014464682127739104*0.0001+z*q2;
      q2 = 2.89247864745380683936*0.000001+z*q2;
      q2 = 6.79019408009981274425*0.000000001+z*q2;
      x1 = z*p2/q2;
    }
  x = x0-x1;
  if( code!=0 )
    {
      x = -x;
    }
  result = x;
  return result;
}
