/*
  A halfspace (infinite) region in space, defined by a plane.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>

#include "dbc.hh"
#include "errors.hh"
#include "Point.hh"

#include "RegionCylinder.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// RegionCylinder:: constructors / destructor
// -------------------------------------------------------------
RegionCylinder::RegionCylinder(const std::string& name,
                               const int id,
                               const Point& axis,
                               const Point& p,
                               double radius,
                               const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::CYLINDER, p.dim(), p.dim(), lifecycle),
    p_(p),
    axis_(axis / norm(axis)),
    rad2_(radius * radius)
{
  if (p_.dim() != axis_.dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in point and axis dimensions of RegionCylinder " << Region::get_name();
    Exceptions::amanzi_throw(mesg);
  }
}


// -------------------------------------------------------------
// RegionCylinder::inside -- check if point is inside cylinder
// -------------------------------------------------------------
bool
RegionCylinder::inside(const Point& p) const
{
  auto dp = p - p_;
  double d = L22(dp);

  double prj = axis_ * dp;
  double dist2 = d - prj * prj;

  return dist2 > rad2_ + tol_ ? false : true;
}

} // namespace AmanziGeometry
} // namespace Amanzi
