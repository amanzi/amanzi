/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   BoxRegion.hh
 * @author Rao Garimella, William A. Perkins
 * @date Wed Sep 28 08:54:19 2011
 * 
 * @brief  Declaration of BoxRegion class (adapted from RectangularRegion)
 * 
 * 
 */

#ifndef _BoxRegion_hh_
#define _BoxRegion_hh_

#include "Region.hh"

namespace Amanzi {
  namespace AmanziGeometry {

// -------------------------------------------------------------
//  class BoxRegion
// -------------------------------------------------------------
/// A rectangular region in space, defined by two points

class BoxRegion : public Region {
public:

  /// Default constructor uses two corner points (order not important).

  BoxRegion(const std::string name, const unsigned int id, const Point& p0, const Point& p1);
  BoxRegion(const char *name, const unsigned int id, const Point& p0, const Point& p1);

  /// Protected copy constructor to avoid unwanted copies.
  BoxRegion(const BoxRegion& old);

  /// Destructor
  ~BoxRegion(void);


  // Type of the region
  inline
  RegionType type() const { return BOX; }

  /// Get the first point defining the region
  inline
  const Point& point0(void) const { return p0_; }

  /// Get the second point defining the region
  inline
  const Point& point1(void) const { return p1_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  /// corners
  inline
  void corners(Point *lo_corner, Point *hi_corner) const
  {
    assert (lo_corner != NULL);
    assert (hi_corner != NULL);

    lo_corner->init(p0_.dim());
    lo_corner->set(p0_);
    hi_corner->init(p1_.dim());
    hi_corner->set(p1_);
  }

  // Is the box degenerate - zero length in one direction?
  bool is_degenerate() const;

protected:
  
  const Point p0_;              /**< one corner of the region  */
  const Point p1_;              /**< the other corner of the region */

  /// Is the specified value between the two values (inclusive, order not important)
  static bool between_(const double& x, const double& x0, const double& x1);

};

/// A smart pointer to BoxRegion instances
//
// typedef Teuchos::RCP<BoxRegion> BoxRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef BoxRegion *BoxRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
