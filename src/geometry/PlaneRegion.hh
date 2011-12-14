/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   PlaneRegion.hh
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Declaration of PlaneRegion class
 * 
 * 
 */

#ifndef _PlaneRegion_hh_
#define _PlaneRegion_hh_

#include "Region.hh"

namespace Amanzi {
  namespace AmanziGeometry {

// -------------------------------------------------------------
//  class PlaneRegion
// -------------------------------------------------------------
/// A planar (infinite) region in space, defined by a point and a plane

class PlaneRegion : public Region {
public:

  /// Default constructor uses point and normal

  PlaneRegion(const std::string name, const unsigned int id, const Point& p, const Point& normal);
  PlaneRegion(const char *name, const unsigned int id, const Point& p, const Point& normal);

  /// Protected copy constructor to avoid unwanted copies.
  PlaneRegion(const PlaneRegion& old);

  /// Destructor
  ~PlaneRegion(void);

  // Type of the region
  inline RegionType type() const { return PLANE; }

  /// Get the point defining the plane
  const Point& point(void) const { return p_; }

  /// Get the normal point defining the plane
  const Point& normal(void) const { return n_; }

  /// Is the specified point inside this region - in this case it
  /// means on the plane

  bool inside(const Point& p) const;

protected:
  
  const Point p_;              /* point on the plane */
  const Point n_;              /* normal to the plane */

};

/// A smart pointer to PlaneRegion instances
// typedef Teuchos::RCP<PlaneRegion> PlaneRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef PlaneRegion *PlaneRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
