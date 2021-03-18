//-*-mode:c++; mode:font-lock;-*-

#include <iomanip>
#include <sstream>

#include "Util.hpp"
#include "Catalog.hpp"

void Catalog::Filter::print(std::ostream& str) const
{
  if(glat_min.first)
    str << "- Galactic Lat         >= " << r2d(glat_min.second) << '\n';
  if(ts_min.first) 
    str << "- TS                   >= " << ts_min.second << '\n';
  if(nassoc_min.first) 
    str << "- Num Assoc            >= " << nassoc_min.second << '\n';
  if(min_neighbour_dist_min.first)
    str << "- Min Neighbour Dist   >= " 
	<< r2d(min_neighbour_dist_min.second) << '\n';
  if(max_inter_assoc_dist_max.first)
    str << "- Max Inter Assoc Dist <= " 
	<< r2d(max_inter_assoc_dist_max.second) << '\n';
  if(max_assoc_prob_min.first)
    str << "- Max Assoc Prob       >= " << max_assoc_prob_min.second << '\n';
  for(std::map<unsigned,double>::const_iterator ib=band_ts_min.begin();
      ib!=band_ts_min.end(); ib++)
    str << "- Band TS[" << ib->first << "]           >= " << ib->second << '\n';
  for(std::map<unsigned,double>::const_iterator ib=band_npred_min.begin();
      ib!=band_npred_min.end(); ib++)
    str << "- Band NPpred[" << ib->first << "]       >= " 
	<< ib->second << '\n';
}

bool Catalog::fillFromFITS(tip::Table::ConstRecord & datum, 
			   const FITSHeader& header,
			   const FITSFillOptions& opt)
{
  datum["nickname"].get(nickname);
  datum["ra"].get(ra);
  datum["dec"].get(dec);
  datum["glon"].get(glon);
  datum["glat"].get(glat);
  datum["test_statistic"].get(ts);
  datum["id_number"].get(nassoc);

  ra   = d2r(ra);
  dec  = d2r(dec);
  glat = d2r(glat);
  glon = d2r(glon);

  std::vector<std::string> id_name;
  std::vector<double> id_prob;
  std::vector<double> id_ra;
  std::vector<double> id_dec;

  datum["id_name"].get(id_name);
  datum["id_probability"].get(id_prob);
  datum["id_ra"].get(id_ra);
  datum["id_dec"].get(id_dec);

  min_neighbour_dist    = 0;
  max_inter_assoc_dist  = 0;
  max_assoc_prob        = 0;
  max_assoc_prob_iassoc = 0;

  assoc.resize(nassoc);
  for(unsigned iassoc=0;iassoc<nassoc;iassoc++)
    {
      assoc[iassoc].name = id_name[iassoc];
      assoc[iassoc].prob = id_prob[iassoc];
      assoc[iassoc].ra   = d2r(id_ra[iassoc]);
      assoc[iassoc].dec  = d2r(id_dec[iassoc]);

      if(id_prob[iassoc]>max_assoc_prob)
	{
	  max_assoc_prob = id_prob[iassoc];
	  max_assoc_prob_iassoc = iassoc;
	}

      for(unsigned jassoc=0;jassoc<iassoc;jassoc++)
	{
	  double d = sphere_dist(assoc[iassoc].ra, assoc[iassoc].dec,
				 assoc[jassoc].ra, assoc[jassoc].dec);
	  if(d>max_inter_assoc_dist)max_inter_assoc_dist = d;
	}
    }

  try
    {
      datum["TS_Bands"].get(band_ts);
    }
  catch(const tip::TipException& o)
    {
    }

  try
    {
      datum["Npred_Bands"].get(band_npred);
    }
  catch(const tip::TipException& o)
    {
    }
  
  return true;
}

#define CUT(testvar,OP,cutvar,cf)					\
  if((cutvar).first && (testvar) OP (cutvar).second)			\
    { if(cut_failed)*cut_failed=(cf); return false; }

bool Catalog::passFilter(const Filter& cuts, Filter::Cut* cut_failed)
{
  CUT(fabs(glat),<,cuts.glat_min,
      Filter::FC_ABS_GLAT);
  CUT(ts,<,cuts.ts_min,
      Filter::FC_TS);
  CUT(nassoc,<,cuts.nassoc_min,
      Filter::FC_NASSOC);
  CUT(min_neighbour_dist,<,cuts.min_neighbour_dist_min,
      Filter::FC_MIN_NEIGH_DIST);
  CUT(max_inter_assoc_dist,>,cuts.max_inter_assoc_dist_max,
      Filter::FC_MAX_ASSOC_DIST);
  CUT(max_assoc_prob,<,cuts.max_assoc_prob_min,
      Filter::FC_MAX_ASSOC_PROB);
  for(std::map<unsigned,double>::const_iterator ib = cuts.band_ts_min.begin();
      ib != cuts.band_ts_min.end(); ib++)
    if(ib->first<band_ts.size() && band_ts[ib->first]<ib->second)
      { if(cut_failed)*cut_failed=Filter::FC_BAND_TS; return false; }
  for(std::map<unsigned,double>::const_iterator ib=cuts.band_npred_min.begin();
      ib != cuts.band_npred_min.end(); ib++)
    if(ib->first<band_npred.size() && band_npred[ib->first]<ib->second)
      { if(cut_failed)*cut_failed=Filter::FC_BAND_NPRED; return false; }
  return true;
}

void Catalog::updateNeighbourDist(std::vector<Catalog>& catalog_vec, 
				  double neighbour_ts)
{
  unsigned ncat = catalog_vec.size();
  for(unsigned icat=0;icat<ncat;icat++)
    {
      catalog_vec[icat].min_neighbour_dist = M_PI;
      for(unsigned jcat=0;jcat<ncat;jcat++)
	if(icat != jcat && catalog_vec[jcat].ts >= neighbour_ts)
	  {
	    double d = sphere_dist(catalog_vec[icat].ra, catalog_vec[icat].dec,
				   catalog_vec[jcat].ra, catalog_vec[jcat].dec);
	    if(d<catalog_vec[icat].min_neighbour_dist)
	      catalog_vec[icat].min_neighbour_dist=d;
	  }
    }
}

void Catalog::print(std::ostream& str, bool print_assoc)
{
  std::ostringstream line;
  line 
    << std::fixed 
    << std::setw(17) << std::left << nickname << ' ' << std::right
    << std::setw(6) << std::setprecision(2) << r2d(ra) << ' '
    << std::setw(6) << std::setprecision(2) << r2d(dec) << ' '
    << std::setw(6) << std::setprecision(2) << r2d(glon) << ' '
    << std::setw(6) << std::setprecision(2) << r2d(glat) << ' '
    << std::setw(9) << std::setprecision(2) << ts << ' '
    << std::setw(1) << nassoc << ' '
    << std::setw(5) << std::setprecision(3) << r2d(min_neighbour_dist) << ' '
    << std::setw(5) << std::setprecision(3) << r2d(max_inter_assoc_dist) << ' '
    << std::setw(6) << std::setprecision(4) << max_assoc_prob << ' ' 
    << std::setw(1) << max_assoc_prob_iassoc << '\n';

  str << line.str();

  if(print_assoc)
    for(unsigned iassoc=0;iassoc<nassoc;iassoc++)
      {
	line.str("");
	line 
	  << "---> " 
	  << std::setw(17) << std::left << assoc[iassoc].name << ' '
	  << std::setw(6) << std::setprecision(4) << assoc[iassoc].prob << ' '
     << std::setw(6) << std::setprecision(2) << r2d(assoc[iassoc].ra) << ' '
     << std::setw(6) << std::setprecision(2) << r2d(assoc[iassoc].dec) << '\n';
	str << line.str();
      }
}

void Catalog::printCatalog(std::vector<Catalog>& catalog_vec,
			   std::ostream& str, bool print_assoc)
{
  unsigned ncat = catalog_vec.size();
  unsigned logncat = 0;
  for(unsigned icat=ncat;icat;icat/=10)logncat++;
  for(unsigned icat=0;icat<ncat;icat++)
    {
      std::ostringstream line;
      line << std::setw(logncat) << std::left << icat+1 << ' ';
      str << line.str();
      catalog_vec[icat].print(str,print_assoc);
    }
}

unsigned Catalog::loadFromFITS(std::vector<Catalog>& catalog_vec,
			       const std::string& filename, 
			       const std::string& tablename,
			       double neighbour_ts,
			       const FITSFillOptions& opt)
{
  unsigned ncat =
    fillVectorFromFITS<Catalog>(catalog_vec, filename, tablename, opt);
  updateNeighbourDist(catalog_vec, neighbour_ts);
  return ncat;
}

unsigned Catalog::filterCatalog(std::vector<Catalog>& catalog_vec,
				const Filter& cuts)
{
  unsigned ncat = catalog_vec.size();
  unsigned icat = 0;
  for(unsigned jcat=0;jcat<ncat;jcat++)
    {
      if(catalog_vec[jcat].passFilter(cuts))
	{
	  if(icat != jcat)catalog_vec[icat]=catalog_vec[jcat];
	  icat++;
	}
    }
  catalog_vec.resize(icat);
  return icat;
}
