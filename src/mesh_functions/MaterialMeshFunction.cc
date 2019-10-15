/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "errors.hh"

#include "MaterialMeshFunction.hh"

namespace Amanzi {
namespace Functions {

/* ******************************************************************
 * Ensure uniqueness of the spec and create the set of IDs contained
 * in the Domain of the spec togher with volume fractions.
 ****************************************************************** */
void
MaterialMeshFunction::AddSpec(const Teuchos::RCP<Spec>& spec)
{
  Teuchos::RCP<Domain> domain = spec->first;
  AmanziMesh::Entity_kind kind = domain->second;
  std::map<AmanziMesh::Entity_ID, double>::iterator it;

  // Loop over regions in the spec, getting their ids and adding to the set.
  Teuchos::RCP<MaterialMesh> mat_mesh = Teuchos::rcp(new MaterialMesh());
  for (RegionList::const_iterator region = domain->first.begin();
       region != domain->first.end();
       ++region) {
    // Get the ids from the mesh by region name and entity kind.
    if (mesh_->valid_set_name(*region, kind)) {
      Kokkos::View<AmanziMesh::Entity_ID*> ids;
      Kokkos::View<double*> vofs;
      mesh_->get_set_entities_and_vofs(
        *region, kind, AmanziMesh::Parallel_type::ALL, ids, &vofs);
      // populating default volume fractions (move this to mesh framework?)
      if (vofs.extent(0) == 0) {
        Kokkos::resize(vofs, ids.extent(0));
        for (int i = 0; i < ids.extent(0); ++i) { vofs(i) = 1.0; }
      }

      for (int i = 0; i < ids.extent(0); ++i) {
        AmanziMesh::Entity_ID id = ids(i);
        it = mat_mesh->find(id);
        if (it == mat_mesh->end()) {
          (*mat_mesh)[id] = vofs(i);
        } else {
          it->second = std::max(it->second + vofs(i), 1.0);
        }
      }
    } else {
      Errors::Message msg;
      msg << "Unknown region in processing mesh function spec: \"" << *region
          << "\", kind=" << kind;
      Exceptions::amanzi_throw(msg);
    }
  }

  // Compare the list of ids in this Spec to all previous Specs to ensure
  // uniqueness.
  Teuchos::RCP<MaterialSpecList> other_specs = material_specs_[kind];
  if (other_specs == Teuchos::null) {
    other_specs = Teuchos::rcp(new MaterialSpecList());
    material_specs_[kind] = other_specs;
  } else {
    for (MaterialSpecList::const_iterator us = other_specs->begin();
         us != other_specs->end();
         ++us) {
      const MaterialMesh& tmp = *(*us)->second;

      for (it = mat_mesh->begin(); it != mat_mesh->end(); ++it) {
        if (tmp.find(it->first) != tmp.end()) {
          Errors::Message msg;
          msg << "Conflicting definitions (overlapping) of id sets.";
          Exceptions::amanzi_throw(msg);
        }
      }
    }
  }
  // If we've gotten this far, things are OK.  Add the spec to the lists.
  spec_list_.push_back(spec);
  other_specs->push_back(Teuchos::rcp(new MaterialSpec(spec, mat_mesh)));
};

} // namespace Functions
} // namespace Amanzi
