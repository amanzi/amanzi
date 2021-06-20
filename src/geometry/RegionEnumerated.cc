/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A region enumerated as a list of IDs.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
#include <locale>
#include <string>

#include "dbc.hh"
#include "errors.hh"

#include "RegionEnumerated.hh"

namespace Amanzi {
namespace AmanziGeometry {

RegionEnumerated::RegionEnumerated(const std::string& name,
                                   const int id,
                                   const std::string& entity_str,
                                   const std::vector<Entity_ID>& ents,
                                   const LifeCycleType lifecycle)
  : Region(name, id, false, ENUMERATED, 0, 0, lifecycle),
    entity_str_(entity_str),
    entities_(ents) {
  std::transform(entity_str_.begin(), entity_str_.end(), entity_str_.begin(),
                 [](unsigned char c) { return std::toupper(c); });
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
