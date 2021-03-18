//-*-mode:c++; mode:font-lock;-*-

#ifndef FT2_HPP
#define FT2_HPP

#include <iostream>

#include "VSAAlgebra.hpp"
#include "FITS.hpp"

class FT2
{
public:
  double t_start;
  double t_stop;
  double zn_ra;
  double zn_dec;
  bool in_saa;
  double scz_ra;
  double scz_dec;
  double scx_ra;
  double scx_dec;
  double pole_ra;
  double pole_dec;
  double livetime;
  unsigned lat_mode;
  unsigned lat_config;
  int data_qual;

  double qso_w;
  double qso_x;
  double qso_y;
  double qso_z;

  class FITSFillOptions
  {
  public:
    FITSFillOptions() { }
  };

  bool fillFromFITS(tip::Table::ConstRecord & datum, const FITSHeader& header,
		    const FITSFillOptions& opt = FITSFillOptions());

  void dumpToOStream(std::ostream& str) const;

  static const char* tableName() { return "sc_data"; }
  static unsigned loadFromFITS(std::vector<FT2>& ft2_vec,
			       const std::string& filename, 
			       const std::string& tablename = tableName(),
			       const FITSFillOptions& opt = FITSFillOptions());

  VERITAS::VSAAlgebra::Vec3D vecScZ() const 
  { return VERITAS::VSAAlgebra::Vec3D::makeLatLon(scz_dec, scz_ra); };
  VERITAS::VSAAlgebra::Vec3D vecScX() const 
  { return VERITAS::VSAAlgebra::Vec3D::makeLatLon(scx_dec, scx_ra); };
  VERITAS::VSAAlgebra::Vec3D vecScY() const 
  { return vecScZ() ^ vecScX(); };
  VERITAS::VSAAlgebra::RotationVec3D vecScRot() const
  { return VERITAS::VSAAlgebra::Vec3D::
      makeRotationQuat(qso_w, qso_x, qso_y, qso_z); }
    
  bool operator<(const FT2& o) const { return t_start<o.t_start; }
};

#endif
