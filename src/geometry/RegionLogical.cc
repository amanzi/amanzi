/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

/*
  A region defined by a logical operation on one or two other regions

  Operations supported on a single region are only the NOT operation

  Operations supported on a pair of regions are UNION, SUBTRACT and
  INTERSECT

*/

#include "dbc.hh"
#include "errors.hh"

#include "RegionLogical.hh"

namespace Amanzi {
namespace AmanziGeometry {

//
// RegionLogical:: constructor
// -------------------------------------------------------------
RegionLogical::RegionLogical(const std::string& name,
                             const int id,
                             const std::string& operation_str,
                             const std::vector<std::string>& component_regions,
                             const LifeCycleType lifecycle)
  : Region(name, id, false, RegionType::LOGICAL, 0, 0, lifecycle),
    operation_(BoolOpType::NOBOOLEAN),
    component_regions_(component_regions)
{
  // Region dimension is set arbitrarily as 0 since the set of
  // entities in the mesh will determine the dimension.
  // 0 should trigger potential errors in the future.

  if (operation_str == "complement")
    operation_ = BoolOpType::COMPLEMENT;
  else if (operation_str == "union")
    operation_ = BoolOpType::UNION;
  else if (operation_str == "intersect")
    operation_ = BoolOpType::INTERSECT;
  else if (operation_str == "subtract")
    operation_ = BoolOpType::SUBTRACT;
  else {
    Errors::Message mesg("Unknown logical operation type requested on regions");
    amanzi_throw(mesg);
  }
}


// -------------------------------------------------------------
// RegionLogical::inside
// -------------------------------------------------------------
bool
RegionLogical::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for logical regions because the check may not "
                       "be implemented for one of its component regions");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
