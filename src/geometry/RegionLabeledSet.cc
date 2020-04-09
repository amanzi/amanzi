/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "dbc.hh"
#include "errors.hh"

#include "RegionLabeledSet.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// RegionLabeledSet:: constructor
// -------------------------------------------------------------
RegionLabeledSet::RegionLabeledSet(const std::string& name, const int id,
                                   const std::string& entity_str,
                                   const std::string& file,
                                   const std::string& format,
                                   const std::string& label,
                                   const LifeCycleType lifecycle)
  : Region(name, id, false, LABELEDSET, 0, 0, lifecycle),
    entity_str_(entity_str),
    file_(file),
    format_(format),
    label_(label)
{}


// -------------------------------------------------------------
// RegionLabeledSet::inside
// -------------------------------------------------------------
bool
RegionLabeledSet::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for labeled sets");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
