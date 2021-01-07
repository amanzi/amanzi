/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! A collection of meshes, indexed from a parent entity set.

/*!

A DomainSet is a collection of meshes, generated from a parent mesh and a set of
entities on that mesh.  For instance, if the parent mesh is a 3D, vertically
structured columnar mesh, one might create a set of meshes, each of which are
1D and correspond to cells of the subsurface mesh.  These column meshes might
be indexed from a set of surface faces (or from the surface mesh's cells).

Alternatively, for e.g. a multi-layer snow model, a DomainSet might be a bunch of
mini meshes constructed at each surface cell.

Note, this does not actually know how to construct theses meshes, it only deals
with their names and their relationship to the parent mesh.  To get/store the
actual mesh, this works with State.

*/
#pragma once

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"

#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

class Mesh;

struct DomainSet {

  // this constructor is for domain sets indexed by an entity + regions
  // (e.g. one for each entity in that collection of regions)
  DomainSet(const Teuchos::RCP<const Mesh>& parent_,
            const std::vector<std::string>& regions_,
            Entity_kind kind_,
            const std::string& name_,
            bool by_region_=true,
            std::vector<std::string>* region_aliases=nullptr);

  using DomainSetMap = std::map<Key,Entity_ID>;

  using const_iterator = DomainSetMap::const_iterator;
  const_iterator begin() const { return meshes.begin(); }
  const_iterator end() const { return meshes.end(); }
  std::size_t size() const { return meshes.size(); }

  std::string get_region(int id) {
    if (by_region) return regions[id];
    else return "";
  }

  Teuchos::RCP<const Mesh> parent;
  std::vector<std::string> regions;
  Entity_kind kind;
  std::string name;
  bool by_region;

  DomainSetMap meshes;
};

} // namespace AmanziMesh
} // namespace Amanzi
