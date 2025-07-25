/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

/*
  A point in space

*/

#include "dbc.hh"
#include "errors.hh"
#include "Point.hh"

#include "RegionPoint.hh"

namespace Amanzi {
namespace AmanziGeometry {

RegionPoint::RegionPoint(const std::string& name,
                         const int id,
                         const Point& p,
                         const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::POINT, 1, p.dim(), lifecycle), p_(p)
{}

// -------------------------------------------------------------
// RegionPoint::inside -- check if input point is coincident with point
// -------------------------------------------------------------
bool
RegionPoint::inside(const Point& p) const
{
  if (p.dim() != p_.dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in dimension of RegionPoint \"" << get_name() << "\" and query point.";
    Exceptions::amanzi_throw(mesg);
  }

  bool result(true);
  for (int i = 0; i != p.dim(); ++i) {
    result = result & (fabs(p[i] - p_[i]) < tol_);
  }
  return result;
}

} // namespace AmanziGeometry
} // namespace Amanzi
