/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   PolygonRegion.hh
 * @author Rao Garimella
 * @date Wed Sep 28 08:54:19 2011
 * 
 * @brief  Declaration of PolygonRegion class 
 * 
 * 
 */

#ifndef _PolygonRegion_hh_
#define _PolygonRegion_hh_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class PolygonRegion
// -------------------------------------------------------------
/// A closed polygonal segment of a plane

class PolygonRegion : public Region {
public:

  /// Default constructor uses two corner points (order not important).

  PolygonRegion(const Set_Name& name, const Set_ID id, 
                const unsigned int npolypoints, 
                const std::vector<Point>& polypoints, 
                const double tolerance=1.0e-8,
                const LifeCycleType lifecycle=PERMANENT,
                const VerboseObject *verbobj=NULL);


  /// Protected copy constructor to avoid unwanted copies.
  PolygonRegion(const PolygonRegion& old);

  /// Destructor
  ~PolygonRegion(void);


  // Type of the region
  inline
  RegionType type() const { return POLYGON; }

  inline 
  unsigned int num_points() const { return num_points_; }

  inline
  const std::vector<Point>& points() const { return points_; }

  inline
  const Point& normal() const { return normal_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:
  
  const unsigned int num_points_;    /* Number of points defining polygon */
  const std::vector<Point> points_;  /* Points of the polygon */
  const double tolerance_;           /* Tolerance for checking containment */
  Point normal_;                     /* Normal to the polygon */
  unsigned int elim_dir_;            /* Coord dir to eliminate while projecting
                                        polygon for in/out tests 
                                        0 - yz, eliminate x coord        
                                        1 - xz, eliminate y coord        
                                        2 - xy, eliminate z coord        */

private:
  void init();
};

/// A smart pointer to PolygonRegion instances
//
// typedef Teuchos::RCP<PolygonRegion> PolygonRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef PolygonRegion *PolygonRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
