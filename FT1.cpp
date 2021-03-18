//-*-mode:c++; mode:font-lock;-*-

#include "FT1.hpp"
#include "Util.hpp"

void FT1::Filter::print(std::ostream& str) const
{
  if(e_min.first)
    str << "- Energy               >= " << e_min.second << '\n';
  if(e_max.first)
    str << "- Energy               <= " << e_max.second << '\n';
  if(z_max.first)
    str << "- Zenith               <= " 
	<< r2d(z_max.second) << '\n';
  if(t_min.first)
    str << "- Theta                >= " << r2d(t_min.second) << '\n';
  if(t_max.first)
    str << "- Theta                <= " << r2d(t_max.second) << '\n';
  if(event_class.first)
    str << "- Event Class          == " 
	<< event_class.second << '\n';
  if(conversion_type.first)
    str << "- Conversion type      == " 
	<< conversion_type.second << '\n';
}

bool FT1::fillFromFITS(tip::Table::ConstRecord & datum, 
		       const FITSHeader& header, 
		       const FITSFillOptions& opt)
{
  datum["energy"].get(energy);
  datum["ra"].get(merit_ra);
  datum["dec"].get(merit_dec);
  datum["theta"].get(theta);
  datum["phi"].get(phi);
  datum["zenith_angle"].get(zenith);
  datum["time"].get(time);
  datum["event_id"].get(event_id);
  datum["run_id"].get(run_id);
  datum["event_class"].get(event_class);
  datum["conversion_type"].get(conversion_type);

  merit_ra  = d2r(merit_ra);
  merit_dec = d2r(merit_dec);
  theta     = d2r(theta);
  phi       = d2r(phi);
  zenith    = d2r(zenith);

  ra        = merit_ra;
  dec       = merit_dec;

  x_dbl.resize(opt.x_dbl_name.size());
  for(unsigned ix=0;ix<opt.x_dbl_name.size();ix++)
    datum[opt.x_dbl_name[ix]].get(x_dbl[ix]);

  return true;
}

#define CUT(testvar,OP,cutvar,cf)					\
  if((cutvar).first && (testvar) OP (cutvar).second)			\
    { if(cut_failed)*cut_failed=(cf); return false; }

bool FT1::passFilter(const Filter& cuts, Filter::Cut* cut_failed)
{
  CUT(energy,<,cuts.e_min,Filter::FC_E_MIN);
  CUT(energy,>,cuts.e_max,Filter::FC_E_MAX);
  CUT(zenith,>,cuts.z_max,Filter::FC_Z_MAX);
  CUT(theta,<,cuts.t_min,Filter::FC_T_MIN);
  CUT(theta,>,cuts.t_max,Filter::FC_T_MAX);
  CUT(event_class,!=,cuts.event_class,Filter::FC_EVENT_CLASS);
  CUT(conversion_type,!=,cuts.conversion_type,Filter::FC_CONVERSION_TYPE);
  return true;
}

unsigned FT1::loadFromFITS(std::vector<FT1>& ft1_vec,
			   const std::string& filename, 
			   const std::string& tablename,
			   const FITSFillOptions& opt)
{
  return fillVectorFromFITS<FT1>(ft1_vec, filename, tablename, opt);  
}

FT1FilterVisitor::~FT1FilterVisitor()
{
  // nothing to see here
}

void FT1FilterVisitor::visitElement(unsigned irow, FT1& t)
{
  if(t.passFilter(m_filter))
    PassThroughFITSVectorVisitor<FT1>::visitElement(irow, t);
}
