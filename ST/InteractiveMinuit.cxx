/** 
 * @file InteractiveMinuit.cxx
 * @brief InteractiveMinuit class implementation
 * @author S. Fegan (based heavily on code by P. Nolan)
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/optimizers/src/Minuit.cxx,v 1.14 2007/11/30 21:21:35 jchiang Exp $
 */

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <limits>
#include <cctype>
#include <cassert>
#include <cfloat>
#include <cmath>

#define IMINUIT_USE_READLINE
#ifdef IMINUIT_USE_READLINE
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include "optimizers/dArg.h"
#include "optimizers/InteractiveMinuit.h"
#include "optimizers/Minuit.h"
#include "optimizers/Parameter.h"
#include "optimizers/Exception.h"
#include "optimizers/OutOfBounds.h"

namespace optimizers 
{

  InteractiveMinuit::InteractiveMinuit(Statistic& stat,
				       const std::string& cmdtxt) 
    : Optimizer(stat), m_commands(), m_nparam(0), 
      m_quality(), m_val(), m_distance(), m_stored_param_values()
  {
    const int i5=5, i6=6, i7=7;
    mninit_(&i5, &i6, &i7);
    
    // Parse command text into individual command strings
    std::string word;
    bool escape = false;
    for(std::string::const_iterator ichar = cmdtxt.begin();
	ichar != cmdtxt.end(); ichar++)
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
	    if(!word.empty())m_commands.push_back(word);
	    word.clear();
	  }
	else
	  {
	    word += *ichar;
	  }
	//std::cout << *ichar << " : " << word << '\n';
      }

    if(!word.empty())m_commands.push_back(word);
  }

  const std::vector<double> & InteractiveMinuit::getUncertainty(bool useBase)
  {
    if (useBase) {
      Optimizer::getUncertainty(useBase);
    }
    return m_uncertainty;
  }

  std::vector<std::string> InteractiveMinuit::getParameterList() const
  {
    std::vector<Parameter> params;
    m_stat->getFreeParams(params);
    std::vector<std::string> pnames(params.size());
    for (unsigned iparam = 0; iparam<params.size(); iparam++)
      pnames[iparam] = params[iparam].getUniqueName();
    return pnames;
  }

  InteractiveMinuit::MinStat InteractiveMinuit::getMinimizationStatus() const
  {
    MinStat ms;
    mnstat_(&ms.fmin, &ms.fedm, &ms.errdef, &ms.npari, &ms.nparx, &ms.istat);
    return ms;
  }

  std::vector<InteractiveMinuit::ValAndErr>
  InteractiveMinuit::getValuesAndErrors() const
  {
    std::vector<Parameter> params;
    m_stat->getFreeParams(params);
    MinStat ms = getMinimizationStatus();
    unsigned nparam = ms.nparx;
    if(nparam != params.size())
      throw Exception("iMINUIT: params.size() != nparam");
    std::vector<InteractiveMinuit::ValAndErr> ve_vec(nparam);
    for(unsigned iparam=0;iparam<nparam;iparam++)
      {
	char buffer[11];
	buffer[10]='\0';
	InteractiveMinuit::ValAndErr& ve(ve_vec[iparam]);
	int j = iparam+1;
	int intvar;
	mnpout_(&j, buffer, &ve.val, &ve.err, &ve.bnd_lo, &ve.bnd_hi, &intvar,
		sizeof(buffer)/sizeof(*buffer)-1);
	ve.uname = params[iparam].getUniqueName();
	ve.pname = params[iparam].getName();
	ve.mname = buffer;
	mnerrs_(&j, &ve.err_plus, &ve.err_minus, &ve.err_parabolic,
		&ve.global_corr);
      }
    return ve_vec;
  }

  int InteractiveMinuit::
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
    mncont_(&fcn, &par1, &par2, &npt_req,
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

  std::vector<std::vector<double> > 
  InteractiveMinuit::getCovarianceMatrix() const
  {
    MinStat ms = getMinimizationStatus();
    int nparam = ms.nparx;
    std::vector<double> elements(nparam*nparam);
    mnemat_(&elements.front(), &nparam);
    std::vector<std::vector<double> > 
      matrix(nparam,std::vector<double>(nparam));
    for(unsigned irow(0), iel(0); irow < nparam; irow++)
      for (unsigned icol(0); icol < nparam; icol++) 
	matrix[irow][icol] = elements[iel++];
    return matrix;
  }

  int InteractiveMinuit::find_min_only(int verbose, double tol, int tolType) 
  {
    return find_min(verbose,tol,tolType);
  }

  int InteractiveMinuit::find_min(int verbose, double tol, int tolType) 
  {
    std::vector<std::string>::const_iterator icommand(m_commands.begin());
    if(icommand != m_commands.end())icommand++;
    
    bool no_init = false;
    bool no_param = false; 
    if(icommand != m_commands.end())
      {
	if(*icommand == ".n")no_init = true;
	else if(*icommand == ".N")no_init = true, no_param = true;
	if(no_init)icommand++;
      }

    typedef std::vector<Parameter>::iterator pptr;

    std::vector<Parameter> params;
    m_stat->getFreeParams(params);
    m_nparam = params.size();
    int errorFlag;

    int minuitVerbose = verbose - 1;
    if (minuitVerbose >= 0) {
      const int i5=5, i6=6, i7=7;
      mintio_(&i5, &i6, &i7);
    }

    if(!no_init)
      {
	std::ostringstream pline;
	pline << "SET PRINT " << minuitVerbose << std::endl;
	doCmd(pline.str()); // Set verbosity of Minuit
	doCmd("SET NOWARN");
	doCmd("SET ERR 0.5");  // Delta value = 1/2: correct for likelihood
	doCmd("SET GRAD 1");  // Use gradient calculated by fcn
      }

    // Tell Minuit about parameter values, names, bounds, etc.
    std::map<std::string,unsigned> pnames;
    for (pptr p = params.begin(); p != params.end(); p++) 
      {
	double scale = 1.0; // Is this the best choice?
	double value = p->getValue();
	double lowerBound = p->getBounds().first;
	double upperBound = p->getBounds().second;
	int j = p - params.begin() + 1;
	if(!no_param)
	  mnparm_(&j, p->getName().c_str(), &value, &scale, 
		  &lowerBound, &upperBound, &errorFlag, p->getName().size());
	pnames[p->getUniqueName()] = j;
      }

    if(verbose)
      std::cout 
	<< std::endl
	<< "iMINUIT: interactive MINUIT optimizer. Type \".h\" for help."
	<< std::endl
	<< std::endl;

    double tolerance = 0.;
    if (tolType == ABSOLUTE)
      tolerance = 2000. * tol;
    else if (tolType == RELATIVE)
      tolerance = 2000. * tol * fabs(m_stat->value());

    // Read in commands and send them to Minuit
    std::istream* input(&std::cin);
    std::list<std::istream*> input_list;
    bool go = true;
    int last_cmd_status = 0;

    while(go)
      {
	std::string command;
	if(icommand != m_commands.end())
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
	  splitCommand(command, pnames, tolerance);

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
	      case 'N':
	      case 'n':
		std::cout << "iMINUIT: commands \".n\" and \".N\" must be given first in parameter list"
			  << std::endl;
		break;
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
		  throw Exception("InteractiveMinuit exception requested", -1);
		else
		  {
		    std::ostringstream cstream;
		    cstream << cmdwords[1];
		    for(unsigned iword=2; iword<ncmdword; iword++)
		      cstream << ' ' << cmdwords[iword];
		    throw Exception(cstream.str(), -1);
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
	      case 'u':
		if(ncmdword != 2)
		  {
		    std::cout
		      << "iMINUIT: error: .u command needs FILENAME"
		      << std::endl;
		    continue;
		  }
		else
		  {
		    std::ifstream istr(cmdwords[1].c_str());
		    if(istr.good())
		      {
			unDumpState(istr);
			std::cout << "iMINUIT: read state from file: "
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
	      case 't':
		if(last_cmd_status != 4)
		  {
		    cmdwords.erase(cmdwords.begin());
		    goto reruncommand;
		  }
	      case 'T':
		if(last_cmd_status == 4)
		  {
		    cmdwords.erase(cmdwords.begin());
		    goto reruncommand;
		  }
	      case 'S':
		{
		  // Get values and errors from Minuit
		  std::vector<ValAndErr> ve = getValuesAndErrors();
		  m_stored_param_values.resize(ve.size());
		  for(unsigned int ive=0; ive<ve.size(); ive++)
		    m_stored_param_values[ive] = ve[ive].val;
		}    
		break;
	      case 'R':
		if(m_stored_param_values.size() == m_nparam)		    
		  for(unsigned int ive=0; ive<m_nparam; ive++)
		    {
		      std::ostringstream cstream;
		      cstream
			<< "SET PAR " << ive+1 << ' '
			<< std::setprecision(12) << std::scientific
			<< m_stored_param_values[ive];
		      doCmd(cstream.str());
		    }
		break;
	      case 'L':
		{
		  double factor=0.1;
		  if(ncmdword==2)
		    std::istringstream(cmdwords[1]) >> factor;
		  // Get values and errors from Minuit
		  std::vector<ValAndErr> ve = getValuesAndErrors();
		  for(unsigned int ive=0; ive<ve.size(); ive++)
		    {
		      double bhi = ve[ive].bnd_hi;
		      double blo = ve[ive].bnd_lo;
		      double val = ve[ive].val;
		      bool logspace = false;

		      if(blo>0 && bhi>blo*10)
			{
			  logspace = true;
			  bhi = std::log(bhi);
			  blo = std::log(blo);
			  val = std::log(val);
			}

		      double window=(bhi-blo)*factor;
		      double setval = val;

		      if(val > bhi-window)
			setval = bhi-window;
		      else if(val < blo+window)
			setval = blo+window;

		      if(setval != val)
			{
			  if(logspace)setval = std::exp(setval);
			  std::ostringstream cstream;
			  cstream
			    << "SET PAR " << ive+1 << ' '
			    << std::setprecision(12) << std::scientific
			    << setval;
			  doCmd(cstream.str());
			}
		    }
		}    
		break;

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

    // Get convergence status
    MinStat ms = getMinimizationStatus();

    // Get values and errors from Minuit
    std::vector<ValAndErr> ve = getValuesAndErrors();

    // Put new parameter values back into the Statistic
    std::vector<double> paramValues(ve.size());
    for(unsigned int ive=0; ive<ve.size(); ive++)
      paramValues[ive] = ve[ive].val;
    m_stat->setFreeParamValues(paramValues);

    // Get parameter uncertainties
    m_uncertainty.resize(ve.size());
    for(unsigned int ive=0; ive<ve.size(); ive++)
      m_uncertainty[ive] = ve[ive].err_parabolic;
    
    m_nparam=0;
  } // End of find_min

  std::pair<double,double> InteractiveMinuit::Minos(unsigned int n) {
    integer npar(m_stat->getNumFreeParams());
    if (n < 1 || n > npar) {
      throw Exception("Parameter number out of range in Minos", n);
    }
    std::ostringstream mcmd;
    mcmd << "MINOS 0 " << n;
    m_nparam = npar;
    doCmd(mcmd.str());
    double eplus, eminus, eparab, globcc;
    integer my_n = n;
    mnerrs_(&my_n, &eplus, &eminus, &eparab, &globcc);
    m_nparam = 0;
    return std::pair<double,double>(eminus,eplus);
  }

  std::vector<std::vector<double> > 
  InteractiveMinuit::covarianceMatrix() const 
  {
    return getCovarianceMatrix();
  }
  
  std::ostream& InteractiveMinuit::put (std::ostream& s) const
    {
      s << "MINUIT returned a value of " << m_val << std::endl;
      s << "and an estimated distance of " << m_distance << std::endl;
      return s;
    }

  // **************************************************************************
  //
  // PRIVATE FUNCTIONS
  //
  // **************************************************************************

  int InteractiveMinuit::doCmd(std::string command) 
  {
    // Pass a command string to Minuit
    int errorFlag = 0;
    void * vthis = static_cast<void *>(this);
    mncomd_(&fcn, command.c_str(), &errorFlag, vthis,
	    command.length());
    return errorFlag;
  }

  void InteractiveMinuit::
  fcn(int* npar, double* grad, double* fcnval,
      double* xval, int* iflag, void* futil)
  {
    // What a hack!  Minuit thinks futil is a function pointer.  It's
    // been hijacked to be a pointer to InteractiveMinuit, so this
    // static function can use it.
    InteractiveMinuit* im = static_cast<InteractiveMinuit *>(futil);

    // This is the function that Minuit minimizes
    std::vector<double> parameters(xval, xval+im->m_nparam);
    
    Statistic* statp = &im->stat();
    statp->setFreeParamValues(parameters);
    
    *fcnval = -statp->value();
    if (*iflag == 2)
      { 
	// Return gradient values
	std::vector<double> gradient;
	statp->getFreeDerivs(gradient);
	for (int i=0; i < *npar; i++)
	  grad[i] = -gradient[i];
      }
  }
  
  std::vector<std::string> InteractiveMinuit::
  splitCommand(const std::string& cmdtxt,
	       const std::map<std::string,unsigned>& parameters,
	       double tolerance) const
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
	    vstream << tolerance;
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

  void InteractiveMinuit::printHelp(bool long_help) const
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
.N     - do not run any of the standard initialization commands, and do\n\
         not specify the list of parameters to MINUIT. This option can \n\
         ONLY be specified on the command line. It is probably not\n\
         useful.\n\
\n\
.n     - do not run any of the standard initialization commands, but DO\n\
         pass the list of parameters to MINUIT. This option can ONLY\n\
         be specified on the command line. It MAY be useful in some\n\
         circumstances.\n\
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
.u XXX - undump the parameter values from a previous dump.\n\
\n\
.c X Y - compute a 2-D error contour using MINOS. X and Y are the variables\n\
         to use on the X and Y axes of the contour respectively.\n\
\n\
.x/.q  - terminate minimization session and allow ScienceTools to\n\
         complete its work (such as calculating TS values in gtlike).\n\
         This should only be done when you are happy you have\n\
         optimized the model and calculated the errors.\n\
\n\
.S/.R  - store current paramater values (\".S\") in permanent storage that\n\
         survives across calls to minimizer, or restore them (\".R\")\n\
\n\
.L XXX - set values that are within a factor XXX (default 0.1) of the limits\n\
         to be outside that range\n\
\n\
.X XXX - throw an exception with message \"XXX\".\n\
\n\
.t XXX - execute the command \"XXX\" only if the last Minuit command ended\n\
.T XXX   in success (\".t\") or in failure (\".T\").\n";
    if(long_help)
      std::cout << "\
\n\
In addition to these commands iMINUIT will replace some words starting\n\
with a \"$\". \"$TOL\" is replaced with the value of the tolerance that\n\
ScienceTools would have used. \"$CHI2(P,N)\" is translated into the\n\
difference in the likelihood function that corresponds to a probability\n\
level of P with N degrees of freedom. Finally, \"$XXX\", where XXX is the\n\
 name of one of the parameters (listed by \".l\") is replaced with the\n\
number of the parameter that MINUIT will use.\n\
\n\
For example, to emulate the standard ScienceTools MINUIT optimizer issue\n\
the following commands:\n\
\n\
iMINUIT> MIGRAD 200 $TOL\n\
iMINUIT> HESSE\n\
iMINUIT> .q\n\
\n\
A more complex example, using MINOS to estimate the errors of two\n\
parameters, and plotting a contour of the 68% region defined by the two\n\
of them is given below:\n\
\n\
iMINUIT> SET STRATEGY 2       # use strictest convergence and error calc\n\
iMINUIT> MIGRAD 0 $TOL        # run MIGRAD with no call limit\n\
iMINUIT> HESSE                # calculate error matrix\n\
iMINUIT> .l                   # get list of parameters\n\
iMINUIT> MINOS 0 4 5          # run MINOS on parameters 4&5, with no call limit\n\
iMINUIT> .D fit.dat           # dump parameter values, errors and corr mat\n\
iMINUIT> SET ERR $CHI2(.68,2) # set error level to 68% with 2 DOF\n\
iMINUIT> SET PRINT 1          # set verbosity so contour points are printed\n\
iMINUIT> MNC 4 5              # calculate confidence contour on parms. 4&5\n\
iMINUIT> .q                   # return to ScienceTools\n\
\n\
If you have worked out a program that you are happy with, and would like\n\
to apply repeatedly, you can enter the program directly on the command line\n\
rather than typing it in interactively. This is done by appending the\n\
commands to the name of the minimizer, each line separated by a comma.\n\
For example, set the minimzer name to:\n\
\n\
   'INTERACTIVEMINUIT,MIGRAD 200 $TOL,HESSE,.q'\n\
\n\
to run the program given in the first example above.\n\
\n\
More information is available on Confluence. Read the CERN MINUIT manual\n\
for lots of advice on how to successfully solve difficult optimization\n\
problems.\
" << std::endl;
  }

  void InteractiveMinuit::printParameters() const
  {
    std::vector<std::string> pnames = getParameterList();
    std::cout << "List of parameters:" << std::endl;
    for (unsigned iparam = 0; iparam<pnames.size(); iparam++)
      std::cout << iparam+1 << ' ' << pnames[iparam] << std::endl;
  }

  //! Read parameter values back from dump of previos state
  void InteractiveMinuit::unDumpState(std::istream& str)
  {
    std::string line;
    std::getline(str, line);
    while(str)
      {
	std::istringstream lstr(line);
	std::string cmd;
	lstr >> cmd;
	if(cmd == "iMINUIT_VE")
	  {
	    unsigned ipar;
	    std::string name;
	    std::string blo;
	    std::string bhi;
	    std::string val;
	    lstr >> ipar >> name >> blo >> bhi >> val;
	    std::ostringstream cstream;
	    cstream
	      << "SET PAR " << ipar << ' ' << val;
	    doCmd(cstream.str());
	  }
	std::getline(str, line);
      }
  }

  //! Print dump of current state
  void InteractiveMinuit::dumpState(std::ostream& str, 
				    bool dumpCovarianceMatrix) const
  {
    MinStat ms = getMinimizationStatus();

    {
      std::ostringstream lstr;
      int np = 3;
      if(ms.fmin>0)np = unsigned(ceil(log10(ms.fmin)))+3;
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
      {
	if(ve[ive].uname.size()>maxusize)maxusize=ve[ive].uname.size();
	if(ve[ive].pname.size()>maxpsize)maxpsize=ve[ive].pname.size();
      }

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
	     << std::setw(maxusize) << std::left << ve[ive].uname << ' '
	  // << std::setw(maxpsize) << std::left << ve[ive].pname << ' ' 
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
	  InteractiveMinuit::getCovarianceMatrix();
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

  double InteractiveMinuit::invchisquaredistribution(double v, double y)
  {
    double result;
    if(y<0 || y>1 || v<1)
      throw Exception("Domain error in InvChiSquareDistribution");
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

  double InteractiveMinuit::invincompletegammac(double a, double y0)
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
        d = (a-1)*log(x)-x-lgm;
        if( d<-709.78271289338399 )
	  {
            d = 0.0625;
            break;
	  }
        d = -exp(d);
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

  double InteractiveMinuit::incompletegammac(double a, double x)
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
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax<-709.78271289338399 )
      {
        result = 0;
        return result;
      }
    ax = exp(ax);
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

  double InteractiveMinuit::incompletegamma(double a, double x)
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
    ax = a*log(x)-x-lngamma(a, tmp);
    if( ax<-709.78271289338399 )
      {
        result = 0;
        return result;
      }
    ax = exp(ax);
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

  double InteractiveMinuit::lngamma(double x, double& sgngam)
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
        z = q*sin(M_PI*z);
        result = logpi-log(z)-w;
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
            result = log(z);
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
        result = log(z)+p;
        return result;
      }
    q = (x-0.5)*log(x)-x+ls2pi;
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

  double InteractiveMinuit::invnormaldistribution(double y0)
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
    x = sqrt(-2.0*log(y));
    x0 = x-log(x)/x;
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

} // namespace optimizers
