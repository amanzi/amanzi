/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "boost/algorithm/string.hpp"

#include "dbc.hh"
#include "errors.hh"

#include "RegionEnumerated.hh"

namespace Amanzi {
namespace AmanziGeometry {

RegionEnumerated::RegionEnumerated(const std::string& name, const int id,
                                   const std::string& entity_str,
                                   const std::vector<Entity_ID>& ents,
                                   const LifeCycleType lifecycle)
  : Region(name, id, false, ENUMERATED, 0, 0, lifecycle),
    entity_str_(boost::to_upper_copy<std::string>(entity_str)),
    entities_(ents)
{
  // Region dimension is set arbitrarily as 0 since the set of
  // entities in the mesh will determine the dimension
}


// -------------------------------------------------------------
// EnumeratedSetRegion::inside
// -------------------------------------------------------------
bool
RegionEnumerated::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for enumerated sets");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
