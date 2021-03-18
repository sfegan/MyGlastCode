//-*-mode:c++; mode:font-lock;-*-

#include "IRF.hpp"
#include "VSFileUtility.hpp"

using namespace VERITAS;

// ****************************************************************************
//
// CALIBRATION DATABASE
//
// ****************************************************************************

bool IRFCalDB::fillFromFITS(tip::Table::ConstRecord & datum,
			    const FITSHeader& header,
			    const FITSFillOptions& opt)
{
  datum["telescop"].get(telescope);
  datum["instrume"].get(instrument);
  datum["detnam"].get(detector);
  datum["filter"].get(filter);
  datum["cal_dev"].get(cal_dev);
  datum["cal_dir"].get(cal_directory);
  datum["cal_file"].get(cal_file);
  datum["cal_clas"].get(cal_class);
  datum["cal_dtyp"].get(cal_dataset_type);
  datum["cal_cnam"].get(cal_code_name);
  datum["cal_cbd"].get(cal_param_limits);
  datum["cal_xno"].get(cal_xno);
  datum["cal_vsd"].get(cal_start_date);
  datum["cal_vst"].get(cal_start_time);
  datum["ref_time"].get(ref_time);
  datum["cal_qual"].get(cal_quality);
  datum["cal_date"].get(cal_date);
  datum["cal_desc"].get(cal_description);

  return true;
}

unsigned IRFCalDB::loadFromFITS(std::vector<IRFCalDB>& caldb_vec,
				const std::string& filename, 
				const std::string& tablename,
				const FITSFillOptions& opt)
{
  return fillVectorFromFITS(caldb_vec, filename, tablename, opt);  
}

std::string IRFCalDB::defaultCalDB()
{
  return "data/glast/lat/caldb.indx";
}

// ****************************************************************************
//
// BASE FOR 2D TABLES
//
// ****************************************************************************

bool IRFTable2D::fillFromFITS(tip::Table::ConstRecord & datum, 
			      const FITSHeader& header,
			      const FITSFillOptions& opt)
{
  datum["energ_lo"].get(elo);
  datum["energ_hi"].get(ehi);
  datum["ctheta_lo"].get(ctlo);
  datum["ctheta_hi"].get(cthi);

  xc.resize(elo.size());
  for(unsigned ie=0;ie<elo.size();ie++)
    xc[ie]=std::log10(std::sqrt(elo[ie]*ehi[ie]));
  
  yc.resize(ctlo.size());
  for(unsigned ict=0;ict<ctlo.size();ict++)yc[ict]=0.5*(ctlo[ict]+cthi[ict]);

  return (elo.size() == ehi.size())&&(ctlo.size() == cthi.size());
}

// ****************************************************************************
//
// EFFECTIVE AREA
//
// ****************************************************************************


bool IRFEffArea::fillFromFITS(tip::Table::ConstRecord & datum, 
			      const FITSHeader& header,
			      const FITSFillOptions& opt)
{
  bool basefill = IRFTable2D::fillFromFITS(datum, header);
  datum["effarea"].get(effarea);
  return basefill && (effarea.size() == elo.size()*ctlo.size());
}

unsigned IRFEffArea::loadFromFITS(std::vector<IRFEffArea>& ea_vec,
				  const std::string& filename, 
				  const std::string& tablename,
				  const FITSFillOptions& opt)
{
  return fillVectorFromFITS(ea_vec, filename, tablename, opt);
}

// ****************************************************************************
//
// EFFICIENCY FACTOR
//
// ****************************************************************************

bool IRFEfficiencyPars::fillFromFITS(tip::Table::ConstRecord & datum, 
				     const FITSHeader& header,
				     const FITSFillOptions& opt)
{
  std::vector<double> pars;
  datum["efficiency_pars"].get(pars);
  if(pars.size() != 6)return false;
  a0        = pars[0];
  b0        = pars[1];
  a1        = pars[2];
  log10e_b1 = pars[3];
  a2        = pars[4];
  log10e_b2 = pars[5];
  b1 = (a0-a1)*log10e_b1+b0;
  b2 = (a1-a2)*log10e_b2+b1;
  return true;
}

unsigned IRFEfficiencyPars::
loadFromFITS(std::vector<IRFEfficiencyPars>& ea_vec,
	     const std::string& filename, const std::string& tablename,
	     const FITSFillOptions& opt)
{
  return fillVectorFromFITS(ea_vec, filename, tablename, opt);
}

// ****************************************************************************
//
// PSF
//
// ****************************************************************************

bool IRFPSFScalingParams::fillFromFITS(tip::Table::ConstRecord & datum, 
				       const FITSHeader& header,
				       const FITSFillOptions& opt)
{
  datum["psfscale"].get(c);
  return true;
}

unsigned IRFPSFScalingParams::
loadFromFITS(std::vector<IRFPSFScalingParams>& param_vec,
	     const std::string& filename, const std::string& tablename,
	     const FITSFillOptions& opt)
{
  return fillVectorFromFITS(param_vec, filename, tablename, opt);  
}

bool IRFPSF::fillFromFITS(tip::Table::ConstRecord & datum, 
			  const FITSHeader& header,
			  const FITSFillOptions& opt)
{
  bool basefill = IRFTable2D::fillFromFITS(datum, header);
  datum["gcore"].get(gcore);
  datum["gtail"].get(gtail);
  datum["ncore"].get(ncore);
  try
    {
      // New function - P6_V8 onwards
      datum["ntail"].get(ntail);
      datum["score"].get(score);
      datum["stail"].get(stail);
    }
  catch(const tip::TipException& o)
    {
      datum["sigma"].get(score);
      stail = score;
      ntail.resize(gcore.size());
      for(unsigned i=0;i<gcore.size();i++)
	{
	  const double gc = gcore[i];
	  const double gt = gtail[i];
	  ntail[i] = ((1.0-1.0/gc)*std::pow(1.0+10.0/gc,-gc))
	    /((1.0-1.0/gt)*std::pow(1.0+10.0/gt,-gt));
	}
    }
  return basefill && (gcore.size() == elo.size()*ctlo.size());
}

unsigned IRFPSF::
loadFromFITS(std::vector<IRFPSF>& psf_vec,
	     const std::string& filename, const std::string& tablename,
	     const FITSFillOptions& opt)
{
  return fillVectorFromFITS(psf_vec, filename, tablename, opt);  
}

void IRFPSF::setScalingParams(const IRFPSFScalingParams& p, bool front)
{
  if(front)
    { m_c0=p.c[0]; m_c1=p.c[1]; m_cg=p.c[4]; }
  else
    { m_c0=p.c[2]; m_c1=p.c[3]; m_cg=p.c[4]; } 
}

void IRFPSF::setScalingParams(const std::string& fn, bool front)
{
  std::vector<IRFPSFScalingParams> p;
  IRFPSFScalingParams::loadFromFITS(p, fn);
  setScalingParams(p[0],front);
}

// ****************************************************************************
//
// IRF LOADER
//
// ****************************************************************************

IRFs::~IRFs()
{
  delete m_ea;
  delete m_eff;
  delete m_phi;
  delete m_psf;
}

void IRFs::
getIRFFileNames(std::string& ea_fn, std::string& psf_fn, 
		const std::string& irf, bool front, bool verbose,
		const std::string& _caldb, const std::string& indx_fn)
{
  // --------------------------------------------------------------------------
  // FIRST CHECK THE CUSTOM IRFS LIST
  // --------------------------------------------------------------------------

  const char* custom_irf_dir = getenv("CUSTOM_IRF_DIR");
  const char* custom_irf_names = getenv("CUSTOM_IRF_NAMES");
  if(custom_irf_dir && custom_irf_names)
    {
      std::vector<std::string> names;
      VSDataConverter::fromString(names, custom_irf_names);
      for(std::vector<std::string>::const_iterator iname = names.begin();
	  iname != names.end(); iname++)
	if(*iname == irf)
	  {
	    ea_fn = custom_irf_dir;
	    ea_fn += "/aeff_";
	    ea_fn += irf;
	    if(front)ea_fn += "_front.fits";
	    else ea_fn += "_back.fits";

	    psf_fn = custom_irf_dir;
	    psf_fn += "/aeff_";
	    psf_fn += irf;
	    if(front)psf_fn += "_front.fits";
	    else psf_fn += "_back.fits";

	    return;
	  }
    }

  // --------------------------------------------------------------------------
  // READ IN IRFS CALDB AND FIND NAMES OF IRF FILES
  // --------------------------------------------------------------------------

  std::string caldb = _caldb;
  VSFileUtility::expandFilename(caldb,true);

  std::string caldb_fn(caldb+"/"+indx_fn);
  if(verbose)
    std::cout << "Load: CALDB      - " << caldb_fn << '\n';

  std::vector<IRFCalDB> caldb_vec;
  IRFCalDB::loadFromFITS(caldb_vec, caldb_fn);

  std::string verirf("VERSION(");
  verirf += irf;
  verirf += ')';

  ea_fn = "";
  psf_fn = "";

  for(unsigned i=0;i<caldb_vec.size();i++)
    {
      if((caldb_vec[i].cal_param_limits[0] == verirf)
	 && ((front && caldb_vec[i].detector == "FRONT") // Could use XOR
	     ||(!front && caldb_vec[i].detector != "FRONT")))
	{
	  std::string file = caldb;
	  file += '/';
	  file += caldb_vec[i].cal_directory;
	  file += '/';
	  file += caldb_vec[i].cal_file;

	  if(caldb_vec[i].cal_code_name == "RPSF")
	    psf_fn = file;
	  else if(caldb_vec[i].cal_code_name == "EFF_AREA")
	    ea_fn = file;
	}
    }
}

bool IRFs::hasIRF(const std::string& irf, bool front,
		  const std::string& _caldb, const std::string& indx_fn)
{
  std::string caldb = _caldb;
  VSFileUtility::expandFilename(caldb,true);

  // --------------------------------------------------------------------------
  // READ IN IRFS CALDB AND FIND NAMES OF IRF FILES
  // --------------------------------------------------------------------------

  std::string caldb_fn(caldb+"/"+indx_fn);

  std::vector<IRFCalDB> caldb_vec;
  IRFCalDB::loadFromFITS(caldb_vec, caldb_fn);

  std::string verirf("VERSION(");
  verirf += irf;
  verirf += ')';

  for(unsigned i=0;i<caldb_vec.size();i++)
    {
      if((caldb_vec[i].cal_param_limits[0] == verirf)
	 && ((front && caldb_vec[i].detector == "FRONT") // Could use XOR
	     ||(!front && caldb_vec[i].detector != "FRONT")))
	return true;
    }
  
  return false;
}

bool IRFs::loadEA(const std::string& ea_fn, bool front,
		  bool efficiency_correction, bool verbose)
{
  // --------------------------------------------------------------------------
  // READ IN EFFECTIVE AREA
  // --------------------------------------------------------------------------

  std::vector<IRFEffArea> ea;

  if(verbose)
    std::cout << "Load: IRF (area) - " << ea_fn << '\n';
  IRFEffArea::loadFromFITS(ea, ea_fn);
  m_ea = new IRFEffArea(ea[0]);

  if(efficiency_correction)
    {
      if(verbose)
	std::cout << "Load: IRF (efcy) - " << ea_fn << '\n';
      try
	{
	  m_eff = new IRFEfficiency;
	  m_eff->setPars(ea_fn, front);
	}
      catch(...)
	{
	  delete m_eff;
	  m_eff = 0;
	  if(verbose)
	    std::cout 
	      << "Load: Failed .. livetime correction not available in IRFs\n";
	}
    }
  else
    {
      delete m_eff;
      m_eff = 0;
    }

  return true;
}

bool IRFs::loadPSF(const std::string& psf_fn, bool front, bool verbose)
{
  // --------------------------------------------------------------------------
  // READ IN POINT SPREAD FUNCTION
  // --------------------------------------------------------------------------

  std::vector<IRFPSF> psf;

  if(verbose)
    std::cout << "Load: IRF (psf)  - " << psf_fn << '\n';
  IRFPSF::loadFromFITS(psf, psf_fn);
  psf[0].setScalingParams(psf_fn, front);

  m_psf = new IRFPSF(psf[0]);

  return true;
}

bool IRFs::
loadAllIRFs(const std::string& irf, bool front, bool efficiency_correction, 
	    bool verbose, const std::string& caldb, const std::string& indx_fn)
{
  // GET FILE NAMES FROM THE CALDB
  std::string ea_fn;
  std::string psf_fn;
  getIRFFileNames(ea_fn, psf_fn, irf, front, verbose, caldb, indx_fn);

  bool loaded_all = true;
  // LOAD EA, EFF AND PHI
  if(!ea_fn.empty())loaded_all &= loadEA(ea_fn, front, efficiency_correction, 
					 verbose);
  // LOAD PSF
  if(!psf_fn.empty())loaded_all &= loadPSF(psf_fn, front, verbose);

  return loaded_all;
}
