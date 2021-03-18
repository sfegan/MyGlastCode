//-*-mode:c++; mode:font-lock;-*-

#ifndef CATALOG_HPP
#define CATALOG_HPP

#include <iostream>

#include "FITS.hpp"

#include "VERITAS/VSAAlgebra.hpp"

class Catalog
{
public:

  class Filter
  {
  public:
    Filter(): 
      glat_min(false,0), ts_min(false,0), nassoc_min(false,0),
      min_neighbour_dist_min(false,0), max_inter_assoc_dist_max(false,0),
      max_assoc_prob_min(false,0), band_ts_min(), band_npred_min() { }

    enum Cut { FC_ABS_GLAT, FC_TS, FC_NASSOC, FC_MIN_NEIGH_DIST,
	       FC_MAX_ASSOC_DIST, FC_MAX_ASSOC_PROB, 
	       FC_BAND_TS, FC_BAND_NPRED };
    std::pair<bool,double> glat_min;
    std::pair<bool,double> ts_min;
    std::pair<bool,unsigned> nassoc_min;
    std::pair<bool,double> min_neighbour_dist_min;
    std::pair<bool,double> max_inter_assoc_dist_max;
    std::pair<bool,double> max_assoc_prob_min;
    std::map<unsigned,double> band_ts_min;
    std::map<unsigned,double> band_npred_min;
    bool haveCuts() const 
    {
      return glat_min.first || ts_min.first || nassoc_min.first 
	|| min_neighbour_dist_min.first || max_inter_assoc_dist_max.first
	|| max_assoc_prob_min.first; 
    }
    void print(std::ostream& str) const;
  };

  class Association
  {
  public:
    std::string name;
    double prob;
    double ra;
    double dec;
  };

  std::string nickname;
  double ra;
  double dec;
  double glon;
  double glat;
  double ts;
  std::vector<double> band_ts;
  std::vector<double> band_npred;
  unsigned nassoc;
  std::vector<Association> assoc;
  double min_neighbour_dist;
  double max_inter_assoc_dist;
  double max_assoc_prob;
  unsigned max_assoc_prob_iassoc;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  bool passFilter(const Filter& cuts, Filter::Cut* cut_failed = 0);
  void print(std::ostream& str, bool print_assoc = false);

  static const char* tableName() { return "lat_point_source_catalog"; }
  static void printCatalog(std::vector<Catalog>& catalog_vec, 
			   std::ostream& str, bool print_assoc = false);
  static void updateNeighbourDist(std::vector<Catalog>& catalog_vec,
				  double neighbour_ts = 0);
  static unsigned loadFromFITS(std::vector<Catalog>& catalog_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       double neighbour_ts = 0,
			       const FITSFillOptions& opt = FITSFillOptions());
  static unsigned filterCatalog(std::vector<Catalog>& catalog_vec,
				const Filter& cuts);
};

#endif
