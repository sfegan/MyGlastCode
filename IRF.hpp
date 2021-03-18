//-*-mode:c++; mode:font-lock;-*-

#ifndef IRF_HPP
#define IRF_HPP

#include <iostream>
#include <algorithm>
#include <vector>

#include "FITS.hpp"

// ****************************************************************************
//
// CALIBRATION DATABASE
//
// ****************************************************************************

class IRFCalDB
{
public:
  std::string telescope;
  std::string instrument;
  std::string detector;
  std::string filter;
  std::string cal_dev;
  std::string cal_directory;
  std::string cal_file;
  std::string cal_class;
  std::string cal_dataset_type;
  std::string cal_code_name;
  std::vector<std::string> cal_param_limits;
  int cal_xno;
  std::string cal_start_date;
  std::string cal_start_time;
  double ref_time;
  int cal_quality;
  std::string cal_date;
  std::string cal_description;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  static const char* tableName() { return "cif"; }
  static unsigned loadFromFITS(std::vector<IRFCalDB>& caldb_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());

  static std::string defaultCalDB();
};

// ****************************************************************************
//
// BASE FOR 2D TABLES
//
// ****************************************************************************

class IRFTable2D
{
public:
  std::vector<double> elo;
  std::vector<double> ehi;
  std::vector<double> ctlo;
  std::vector<double> cthi;

  std::vector<double> xc;
  std::vector<double> yc;
  
  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  void xyindex(const double log10e, const double ctheta,
	       unsigned& ix, double& dx, unsigned& iy, double& dy) const
  {
    const unsigned nx = xc.size()-1;
    if(log10e<xc[0])ix=0,dx=0;
    else if(log10e>=xc[nx])ix=nx,dx=0;
    else
      {
	ix = 0;
	if(nx>1) // Try to guess x-index but use upper_bound if we guess wrong
	  {
	    const unsigned nx2=nx>>1;
	    ix = nx2+lrint((log10e-xc[nx2])/(xc[nx2+1]-xc[nx2])-0.5);
	    if(ix>nx)ix=nx;
	  }

	if(log10e<xc[ix] || log10e>=xc[ix+1])
	  {
	    std::vector<double>::const_iterator iix =
	      std::upper_bound(xc.begin(),xc.end(),log10e);
	    ix = iix-xc.begin()-1;
	  }

	assert(log10e>=xc[ix] && log10e<xc[ix+1]);

	dx = (log10e-xc[ix])/(xc[ix+1]-xc[ix]);
      }

    const unsigned ny = yc.size()-1;
    if(ctheta<yc[0])iy=0,dy=0;
    else if(ctheta>=yc[ny])iy=ny,dy=0;
    else
      {
	iy = lrint((ctheta-yc[0])/(yc[ny]-yc[0])*double(ny)-0.5);
	if(ctheta<yc[iy] || ctheta>=yc[iy+1])
	  {
	    std::vector<double>::const_iterator iiy =
	      std::upper_bound(yc.begin(),yc.end(),ctheta);
	    iy = iiy-yc.begin()-1;
	  }

	assert(ctheta>=yc[iy] && ctheta<yc[iy+1]);

	dy = (ctheta-yc[iy])/(yc[iy+1]-yc[iy]);
      }    
  }

  double value(const double log10e, const double ctheta, 
	       const std::vector<double>& data) const
  {
    unsigned ix;
    double dx;
    unsigned iy;
    double dy;
    xyindex(log10e, ctheta, ix, dx, iy, dy);

    double val = 0;
    if(dx==0)
      {
	if(dy!=0)
	  val = data[I(ix,iy)]*(1.0-dy) + data[I(ix,iy+1)]*dy;
	else 
	  val = data[I(ix,iy)];
      }
    else if(dy==0)
      {
	val = data[I(ix,iy)]*(1.0-dx) + data[I(ix+1,iy)]*dx;
      }
    else
      
      {
	val = (data[I(ix,iy)]*(1.0-dy) + data[I(ix,iy+1)]*dy)*(1.0-dx)
	  + (data[I(ix+1,iy)]*(1.0-dy) + data[I(ix+1,iy+1)]*dy)*dx;
      }
    return val;
  }
private:
  unsigned I(unsigned ix, unsigned iy) const { return iy*xc.size()+ix; }
};

// ****************************************************************************
//
// EFFECTIVE AREA
//
// ****************************************************************************

class IRFEffArea: public IRFTable2D
{
public:
  std::vector<double> effarea;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  static const char* tableName() { return "effective area"; }
  static unsigned loadFromFITS(std::vector<IRFEffArea>& ea_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());

  double value(const double log10e, const double ctheta) const
  {
    return IRFTable2D::value(log10e,ctheta,effarea);
  }
};

// ****************************************************************************
//
// EFFICIENCY FACTOR
//
// ****************************************************************************

class IRFEfficiencyPars
{
public:
  double a0;
  double b0;
  double a1;
  double b1;
  double a2;
  double b2;
  double log10e_b1;
  double log10e_b2;
  std::vector<double> pars;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  static const char* tableName() { return "efficiency_params"; }
  static unsigned loadFromFITS(std::vector<IRFEfficiencyPars>& eff_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());

  double value(double log10e) const
  {
    if(log10e<log10e_b1)
      return a0*log10e + b0;
    else if(log10e<log10e_b2)
      return a1*log10e + b1;
    else 
      return a2*log10e + b2;
  }  
};

class IRFEfficiency
{
public:
  IRFEfficiency(): m_have_pars(false), m_p0(), m_p1() { }

  void setPars(const IRFEfficiencyPars& _p0, const IRFEfficiencyPars& _p1)
  {
    m_have_pars=true;
    m_p0=_p0;
    m_p1=_p1;
  }

  void setPars(const std::vector<IRFEfficiencyPars>& p, bool front)
  {
    m_have_pars=true;
    if(front)m_p0 = p.at(0), m_p1 = p.at(1);
    else m_p0 = p.at(2), m_p1 = p.at(3);
  }

  void setPars(const std::string& filename, bool front,
	       const std::string& tablename = IRFEfficiencyPars::tableName())
  {
    std::vector<IRFEfficiencyPars> p;
    IRFEfficiencyPars::loadFromFITS(p, filename, tablename);
    setPars(p, front);
  }

  double p0(double log10e) const
  { 
    if(!m_have_pars)return 0.0; 
    else return m_p0.value(log10e); 
  }

  double p1(double log10e) const
  {
    if(!m_have_pars)return 1.0;
    else return m_p1.value(log10e);
  }

  double value(double livetime_frac, double p0_value, double p1_value) const
  {
    return p0_value*livetime_frac + p1_value;
  }

  double value(double livetime_frac, double log10e) const
  {
    if(!m_have_pars)return 1.0;
    return value(livetime_frac, m_p0.value(log10e), m_p1.value(log10e));
  }

private:
  bool m_have_pars;
  IRFEfficiencyPars m_p0;
  IRFEfficiencyPars m_p1;
};

// ****************************************************************************
//
// PHI DEPENDENCE
//
// ****************************************************************************

class IRFPhiModulation
{
public:
};

// ****************************************************************************
//
// PSF
//
// ****************************************************************************

class IRFPSFScalingParams
{
public:
  std::vector<double> c;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  static const char* tableName() { return "psf_scaling_params"; }
  static unsigned loadFromFITS(std::vector<IRFPSFScalingParams>& param_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());
};

class IRFPSF: public IRFTable2D
{
public:  
  std::vector<double> ncore; // <-- Irrelavant parameter not used below
  std::vector<double> ntail;
  std::vector<double> score;
  std::vector<double> stail;
  std::vector<double> gcore;
  std::vector<double> gtail;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  static const char* tableName() { return "rpsf"; }
  static unsigned loadFromFITS(std::vector<IRFPSF>& param_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());
  
  void setScalingParams(const IRFPSFScalingParams& p, bool front=true);
  void setScalingParams(const std::string& fn, bool front=true);

  double value(const double distance, 
	       const double log10e, const double ctheta) const
  {
    const double scale = 
      std::sqrt(sqr(m_c0*std::pow(10,m_cg*(log10e-2)))+sqr(m_c1));

    const double sc = IRFTable2D::value(log10e,ctheta,score)*scale;
    const double gc = IRFTable2D::value(log10e,ctheta,gcore);
    const double sc2 = sqr(sc);

    const double nt = IRFTable2D::value(log10e,ctheta,ntail);
    const double st = IRFTable2D::value(log10e,ctheta,stail)*scale;
    const double gt = IRFTable2D::value(log10e,ctheta,gtail);
    const double st2 = sqr(st);

    const double d2 = sqr(distance);
    const double uc = 0.5*d2/sc2;
    const double ut = 0.5*d2/st2;

    const double pc = (1.0-1.0/gc)*std::pow(1.0+uc/gc,-gc);
    const double pt = (1.0-1.0/gt)*std::pow(1.0+ut/gt,-gt);

    return (pc+nt*pt)/(2.0*M_PI)/(sc2+nt*st2);
  }

  double integral(const double distance, 
		  const double log10e, const double ctheta) const
  {
    const double scale = 
      std::sqrt(sqr(m_c0*std::pow(10,m_cg*(log10e-2)))+sqr(m_c1));

    const double sc = IRFTable2D::value(log10e,ctheta,score)*scale;
    const double gc = IRFTable2D::value(log10e,ctheta,gcore);
    const double sc2 = sqr(sc);

    const double nt = IRFTable2D::value(log10e,ctheta,ntail);
    const double st = IRFTable2D::value(log10e,ctheta,stail)*scale;
    const double gt = IRFTable2D::value(log10e,ctheta,gtail);
    const double st2 = sqr(st);

    const double d2 = sqr(distance);
    const double uc = 0.5*d2/sc2;
    const double ut = 0.5*d2/st2;

    const double Ic = 1.0-std::pow(1.0+uc/gc,1-gc);
    const double It = 1.0-std::pow(1.0+ut/gt,1-gt);

    return (sc2*Ic+nt*st2*It)/(sc2+nt*st2);
  }

  double rP(const double P, const double log10e, const double ctheta) const
  {
    const double scale = 
      std::sqrt(sqr(m_c0*std::pow(10,m_cg*(log10e-2)))+sqr(m_c1));

    const double sc = IRFTable2D::value(log10e,ctheta,score)*scale;
    const double gc = IRFTable2D::value(log10e,ctheta,gcore);
    const double sc2 = sqr(sc);

    const double nt = IRFTable2D::value(log10e,ctheta,ntail);
    const double st = IRFTable2D::value(log10e,ctheta,stail)*scale;
    const double gt = IRFTable2D::value(log10e,ctheta,gtail);
    const double st2 = sqr(st);

    double dl = 0;
    double dr = M_PI/2;

    while(std::abs(dl-dr)>1e-6)
      {
	const double dm  = 0.5*(dl+dr);
	const double dm2 = sqr(dm);
	const double uc  = 0.5*dm2/sc2;
	const double ut  = 0.5*dm2/st2;
	const double Ic  = 1.0-std::pow(1.0+uc/gc,1-gc);
	const double It  = 1.0-std::pow(1.0+ut/gt,1-gt);
	const double I = (sc2*Ic + nt*st2*It)/(sc2+nt*st2);
	if(I>P)dr = dm;
	else dl = dm;
      }
    
    return 0.5*(dl+dr);
  }

private:
  static inline double sqr(const double& x) { return x*x; }
  double m_c0;
  double m_c1;
  double m_cg;
};

// ****************************************************************************
//
// IRF LOADER
//
// ****************************************************************************

class IRFs
{
public:
  IRFs(): m_ea(), m_eff(), m_phi(), m_psf() { }
  ~IRFs();

  static void 
  getIRFFileNames(std::string& ea_fn, std::string& psf_fn,
		  const std::string& irf, bool front, bool verbose = false,
		  const std::string& caldb = "$CALDB",
		  const std::string& indx_fn = IRFCalDB::defaultCalDB());

  static bool hasIRF(const std::string& irf, bool front,
		     const std::string& caldb = "$CALDB",
		     const std::string& indx_fn = IRFCalDB::defaultCalDB());

  bool loadEA(const std::string& ea_fn, bool front,
	      bool efficiency_correction = true, bool verbose = false);
  bool loadPSF(const std::string& psf_fn, bool front, bool verbose = false);

  bool loadAllIRFs(const std::string& irf, bool front, 
		   bool efficiency_correction = true, bool verbose = false,
		   const std::string& caldb = "$CALDB",
		   const std::string& indx_fn = IRFCalDB::defaultCalDB());

  const IRFEffArea* ea() const { return m_ea; }
  const IRFEfficiency* eff() const { return m_eff; }
  const IRFPhiModulation* phi() const { return m_phi; }
  const IRFPSF* psf() const { return m_psf; }

private:
  IRFs(const IRFs&);
  IRFs& operator=(const IRFs&);

  IRFEffArea*        m_ea;
  IRFEfficiency*     m_eff;
  IRFPhiModulation*  m_phi;
  IRFPSF*            m_psf;
};

#endif
