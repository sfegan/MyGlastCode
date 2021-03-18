//-*-mode:c++; mode:font-lock;-*-

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "VERITAS/VSOptions.hpp"
#include "VERITAS/VSAAlgebra.hpp"

#include "FT1.hpp"
#include "FT2.hpp"
#include "Util.hpp"

// ----------------------------------------------------------------------------
//
// SOURCE FINDER
//
// ----------------------------------------------------------------------------

class FoundSourceVisitor: public FITSVectorVisitor<FT1>
{
public:
  FoundSourceVisitor(): FITSVectorVisitor<FT1>() { }
  virtual ~FoundSourceVisitor();
  virtual void visitFSEvent(unsigned irow, FT1& event,
			    unsigned isrc, double ra, double dec, double dist);
};

class SourceFinder: public PassThroughFITSVectorVisitor<FT1>
{
public:
  SourceFinder(FoundSourceVisitor* visitor, double d_max,
	       const std::vector<double>& cat_ra,
	       const std::vector<double>& cat_dec);
  virtual ~SourceFinder();
  virtual void visitElement(unsigned irow, FT1& event);

  unsigned nFSEventsDispatched() const { return m_nfsevents; }
private:
  FoundSourceVisitor* m_visitor;
  double m_d_max;
  unsigned m_ncat;
  std::vector<double> m_cat_ra;
  std::vector<double> m_cat_dec;
  std::vector<double> m_cat_sra;
  std::vector<double> m_cat_cra;
  std::vector<double> m_cat_sdec;
  std::vector<double> m_cat_cdec;
  unsigned m_nfsevents;
};

// ----------------------------------------------------------------------------
//
// RADEC ORIENTATION
//
// ----------------------------------------------------------------------------

class RADECRotation: public PassThroughFITSVectorVisitor<FT1>
{
public:
  RADECRotation(FITSVectorVisitor<FT1>* visitor, 
		const VERITAS::VSAAlgebra::Vec3D& radec_rot =
		VERITAS::VSAAlgebra::Vec3D())
    : PassThroughFITSVectorVisitor<FT1>(visitor), m_radec_rot(radec_rot) { }
  virtual ~RADECRotation();
  virtual void visitElement(unsigned irow, FT1& event);
private:
  VERITAS::VSAAlgebra::Vec3D m_radec_rot;
};

// ----------------------------------------------------------------------------
//
// SPACECRAFT ORIENTATION
//
// ----------------------------------------------------------------------------

class SOEventVisitor: public FITSVectorVisitor<FT1>
{
public:
  SOEventVisitor(): FITSVectorVisitor<FT1>() { }
  virtual ~SOEventVisitor();
  virtual void visitSOEvent(unsigned irow, FT1& event, 
			    VERITAS::VSAAlgebra::Vec3D& so);
};

class EventRADECRecalc: public SOEventVisitor
{
public:
  EventRADECRecalc(FITSVectorVisitor<FT1>* visitor, 
		   const VERITAS::VSAAlgebra::Vec3D& sc_rot =
		   VERITAS::VSAAlgebra::Vec3D(),
		   std::vector<double> fourier_dtheta = 
		   std::vector<double>())
    : SOEventVisitor(), m_visitor(visitor), m_sc_rot(sc_rot),
      m_fourier_dtheta(fourier_dtheta) { }
  virtual ~EventRADECRecalc();
  virtual void visitSOEvent(unsigned irow, FT1& event, 
			    VERITAS::VSAAlgebra::Vec3D& so);
  virtual void visitFileVector(const std::string& filename,
			       const std::string& tablename,
			       unsigned nrow, FITSHeader& header);
  virtual void leaveFileVector();
private:
  FITSVectorVisitor<FT1>* m_visitor;
  VERITAS::VSAAlgebra::Vec3D m_sc_rot;
  std::vector<double> m_fourier_dtheta;
};

class EODatEntryMatcher: public PassThroughFITSVectorVisitor<FT1>
{
public:
  EODatEntryMatcher(SOEventVisitor* visitor, const std::string& eodat_fn);
  virtual void visitElement(unsigned irow, FT1& event);
  unsigned nMatchedEventsDispatched() const { return m_nmatched; }
private:
  SOEventVisitor* m_visitor;
  std::map<std::pair<unsigned,unsigned>, VERITAS::VSAAlgebra::Vec3D> m_eo_map;
  unsigned m_nmatched;
};

class QuadInterpSOCalculator: 
  public PassThroughFITSVectorVisitor<FT1>
{
public:
  QuadInterpSOCalculator(SOEventVisitor* visitor, const std::vector<FT2>& ft2)
    : PassThroughFITSVectorVisitor<FT1>(visitor),
      m_ft2(ft2), m_ift2(m_ft2.begin()), m_visitor(visitor), m_nsoevents() { }
  virtual ~QuadInterpSOCalculator();
  virtual void visitElement(unsigned irow, FT1& event);
  unsigned nSOEventsDispatched() const { return m_nsoevents; }
private:
  const std::vector<FT2>& m_ft2;
  std::vector<FT2>::const_iterator m_ift2;
  SOEventVisitor* m_visitor;
  unsigned m_nsoevents;
};


#endif
