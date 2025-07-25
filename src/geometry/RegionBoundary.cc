/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov(lipnikov@lanl.gov)
*/

/*
  A region consisting of all entities on the domain boundary.

*/

#include "dbc.hh"
#include "errors.hh"
#include "RegionBoundary.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// Region dimension is set arbitrarily as 0 since the set of
// entities in the mesh will determine the dimension.
// -------------------------------------------------------------
RegionBoundary::RegionBoundary(const std::string& name, const int id, const LifeCycleType lifecycle)
  : Region(name, id, false, RegionType::BOUNDARY, 0, 0, lifecycle) {};


// -------------------------------------------------------------
// Required member function
// -------------------------------------------------------------
bool
RegionBoundary::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for region 'boundary'");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
