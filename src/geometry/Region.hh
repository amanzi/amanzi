/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   Region.hh
 * @author William A. Perkins
 * @date Fri Jul 29 12:40:09 2011
 * 
 * @brief  Declaration of the abstract Region class 
 * 
 * 
 */

#ifndef _Region_hh_
#define _Region_hh_

#include <vector>
#include <Teuchos_RCP.hpp>

#include "Point.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class Region
// -------------------------------------------------------------
/// An class that represent a geometric region
/**
 * A Region is just some arbitrary subset of space, that can be
 * specified in a myriad of ways.  At a minimum, there is a need to be
 * able to determine if a point is inside that space.  Other needs to
 * be added later.
 * 
 */
class Region {
public:

  /// Default constructor.
  Region(void);

  /// Copy constructor 
  Region(const Region& old);

  /// Destructor
  virtual ~Region(void);

  /// Is the the specified point inside the Region
  virtual bool inside(const Point& p) const = 0;

  /// Get the extents of the Region
  virtual void extents(Point *pmin, Point *pmax) const = 0;
};

/// A smart pointer to Region instances
typedef Teuchos::RCP<Region> RegionPtr;

/// A thing to hold some Region
typedef std::vector< RegionPtr > RegionVector;

} // namespace AmanziGeometry
} // namespace Amanzi

#endif

