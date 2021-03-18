//-*-mode:c++; mode:font-lock;-*-

#include "FT2.hpp"
#include "Util.hpp"

#include "VSDataConverter.hpp"

using namespace VERITAS;

bool FT2::fillFromFITS(tip::Table::ConstRecord & datum,
		       const FITSHeader& header,
		       const FITSFillOptions& opt)
{
  datum["start"].get(t_start);
  datum["stop"].get(t_stop);
  datum["ra_zenith"].get(zn_ra);
  datum["dec_zenith"].get(zn_dec);
  datum["in_saa"].get(in_saa);
  datum["ra_scz"].get(scz_ra);
  datum["dec_scz"].get(scz_dec);
  datum["ra_scx"].get(scx_ra);
  datum["dec_scx"].get(scx_dec);
  datum["ra_npole"].get(pole_ra);
  datum["dec_npole"].get(pole_dec);
  datum["livetime"].get(livetime);
  datum["lat_mode"].get(lat_mode);
  datum["lat_config"].get(lat_config);
  datum["data_qual"].get(data_qual);

  datum["qsj_1"].get(qso_x);
  datum["qsj_2"].get(qso_y);
  datum["qsj_3"].get(qso_z);
  datum["qsj_4"].get(qso_w);

  zn_ra     = d2r(zn_ra);
  zn_dec    = d2r(zn_dec);
  scz_ra    = d2r(scz_ra);
  scz_dec   = d2r(scz_dec);
  scx_ra    = d2r(scx_ra);
  scx_dec   = d2r(scx_dec);
  pole_ra   = d2r(pole_ra);
  pole_dec  = d2r(pole_dec);

  //scz = Vec3D::makeLatLong(scz_dec, scz_ra);
  //scx = Vec3D::makeLatLong(scx_dec, scx_ra);
  //r =   Vec3D::makeRotationQuat(qso_w,qso_x,qso_y,qso_z);

  return true;
}

unsigned FT2::loadFromFITS(std::vector<FT2>& ft2_vec,
			   const std::string& filename, 
			   const std::string& tablename,
			   const FITSFillOptions& opt)
{
  return fillVectorFromFITS<FT2>(ft2_vec, filename, tablename, opt);  
}

void FT2::dumpToOStream(std::ostream& str) const
{
  str << VSDataConverter::toString(t_start) << ' ' 
      << VSDataConverter::toString(t_stop) << ' '
      << r2d(zn_ra) << ' ' << r2d(zn_dec) << ' ' << in_saa << ' '
      << r2d(scz_ra) << ' ' << r2d(scz_dec) << ' ' 
      << r2d(scx_ra) << ' ' << r2d(scx_dec) << ' '
      << r2d(pole_ra) << ' ' << r2d(pole_dec) << ' '
      << qso_w << ' ' << qso_x << ' ' << qso_y << ' ' << qso_z << '\n';

}
