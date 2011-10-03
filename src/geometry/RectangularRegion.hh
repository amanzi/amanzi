/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   RectangularRegion.hh
 * @author William A. Perkins
 * @date Wed Sep 28 08:54:19 2011
 * 
 * @brief  Declaration of RectangularRegion class
 * 
 * 
 */

#ifndef _RectangularRegion_hh_
#define _RectangularRegion_hh_

#include "Region.hh"

namespace Amanzi {
  namespace AmanziGeometry {

// -------------------------------------------------------------
//  class RectangularRegion
// -------------------------------------------------------------
/// A rectangular region in space, defined by two points

class RectangularRegion : public Region {
public:

  /// Default constructor uses to corner points (order not important).
  RectangularRegion(const Point& p0, const Point& p1);

  /// Protected copy constructor to avoid unwanted copies.
  RectangularRegion(const RectangularRegion& old);

  /// Destructor
  ~RectangularRegion(void);

  /// Get the first point defining the region
  const Point& point0(void) const { return p0_; }

  /// Get the second point defining the region
  const Point& point1(void) const { return p1_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:
  
  const Point p0_;              /**< one corner of the region  */
  const Point p1_;              /**< the other corner of the region */

  /// Is the specified value between the two values (inclusive, order not important)
  static bool between_(const double& x, const double& x0, const double& x1);

};

/// A smart pointer to RectangularRegion instances
typedef Teuchos::RCP<RectangularRegion> RectangularRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
