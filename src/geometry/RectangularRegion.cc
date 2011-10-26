/**
 * @file   RectangularRegion.cc
 * @author William A. Perkins
 * @date Fri Jul 29 12:28:10 2011
 * 
 * @brief  Implementation of RectangularRegion class
 * 
 * 
 */

#include "RectangularRegion.hh"
#include "dbc.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class RectangularRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// RectangularRegion:: constructors / destructor
// -------------------------------------------------------------
RectangularRegion::RectangularRegion(const Point& p0, const Point& p1)
  : Region(), p0_(p0), p1_(p1)
{
  ASSERT(p0_.dim() == p1_.dim());
}

RectangularRegion::RectangularRegion(const RectangularRegion& old)
  : Region(old), p0_(old.p0_), p1_(old.p1_)
{
  // empty
}

RectangularRegion::~RectangularRegion(void)
{
  
}

// -------------------------------------------------------------
// RectangularRegion::between_
// -------------------------------------------------------------
bool
RectangularRegion::between_(const double& x, const double& x0, const double& x1)
{
  return  (x0 >= x && x >= x1) || (x1 >= x && x >= x0);
}

// -------------------------------------------------------------
// RectangularRegion::inside
// -------------------------------------------------------------
bool
RectangularRegion::inside(const Point& p) const
{
  ASSERT(p.dim() == p0_.dim());
  bool result(true);
  for (int i = 0; i < p.dim(); ++i) {
    result = result && between_(p[i], p0_[i], p1_[i]);
  }
  return result;
}

} // namespace AmanziGeometry
} // namespace Amanzi
