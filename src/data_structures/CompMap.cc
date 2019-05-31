/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

  Implementation for CompMap, a dictionary from keys to their maps.
*/

#include "DataStructuresHelpers.hh"
#include "CompMap.hh"

namespace Amanzi {

bool CompMap::SameAs(const CompMap& other) const {
  if (names_ != other.names_) return false;
  for (const auto& name : names_) {
    if (!Amanzi::SameAs(*ComponentMap(name), *other.ComponentMap(name))) return false;
  }
  return true;
}



} // namespace Amanzi
