//-*-mode:c++; mode:font-lock;-*-

#ifndef FT1_HPP
#define FT1_HPP

#include <iostream>
#include <string>
#include <vector>

#include "FITS.hpp"

class FT1
{
public:
  class Filter
  {
  public:
    enum Cut { FC_E_MIN, FC_E_MAX, FC_Z_MAX, FC_T_MIN, FC_T_MAX,
	       FC_EVENT_CLASS, FC_CONVERSION_TYPE };
    Filter()
      : e_min(), e_max(), z_max(), t_min(), t_max(),
	event_class(), conversion_type() { }
    std::pair<bool,double> e_min;
    std::pair<bool,double> e_max;
    std::pair<bool,double> z_max;
    std::pair<bool,double> t_min;
    std::pair<bool,double> t_max;
    std::pair<bool,int> event_class;
    std::pair<bool,int> conversion_type;
    bool haveCuts() const 
    {
      return e_min.first || e_max.first 
	|| z_max.first || t_min.first || t_max.first 
	|| event_class.first || conversion_type.first;
    }
    void print(std::ostream& str) const;
  };

  double energy;
  double merit_ra;
  double merit_dec;  
  double theta;
  double phi;
  double zenith;
  double time;
  unsigned event_id;
  unsigned run_id;
  int event_class;
  int conversion_type;

  double ra;
  double dec;  

  std::vector<double> x_dbl;

  class FITSFillOptions
  {
  public:
    FITSFillOptions(): x_dbl_name() { }
    FITSFillOptions(const std::string& name): 
      x_dbl_name() { x_dbl_name.push_back(name); }
    FITSFillOptions(const std::vector<std::string>& _x_dbl_name): 
      x_dbl_name(_x_dbl_name) { }
    std::vector<std::string> x_dbl_name;
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());
  bool passFilter(const Filter& cuts, Filter::Cut* cut_failed = 0);

  static const char* tableName() { return "events"; }
  static unsigned loadFromFITS(std::vector<FT1>& ft1_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());
};

class FT1FilterVisitor: public PassThroughFITSVectorVisitor<FT1>
{
public:
  FT1FilterVisitor(FITSVectorVisitor<FT1>* visitor, 
		   const FT1::Filter& filter):
    PassThroughFITSVectorVisitor<FT1>(visitor), m_filter(filter) { }
  virtual ~FT1FilterVisitor();
  virtual void visitElement(unsigned irow, FT1& t);
private:
  FT1::Filter m_filter;
};

#endif
