/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "errors.hh"

#include "UniqueMeshFunction.hh"

namespace Amanzi {
namespace Functions {

// Overload the AddSpec method to check uniqueness.
void
UniqueMeshFunction::AddSpec(const Teuchos::RCP<Spec>& spec)
{
  // Ensure uniqueness of the spec and create the set of IDs contained in
  // the domain of the spec.

  Teuchos::RCP<Domain>& domain = spec->first;
  AmanziMesh::Entity_kind kind = domain->second;

  // Loop over regions in the spec, getting their ids and adding to the set.
  Teuchos::RCP<MeshIDs> this_spec_ids = Teuchos::rcp(new MeshIDs());
  for (auto region = domain->first.begin(); region != domain->first.end();
       ++region) {
    // Get all region IDs by the region name and entity kind.
    if (mesh_->valid_set_name(*region, kind)) {
      Kokkos::View<AmanziMesh::Entity_ID*> id_list;
      mesh_->get_set_entities(
        *region, kind, AmanziMesh::Parallel_type::ALL, id_list);
      for (int i = 0; i < id_list.extent(0); ++i)
        this_spec_ids->insert(id_list(i));
    } else {
      std::stringstream m;
      m << "Unknown region in processing mesh function spec: \"" << *region
        << "\"";
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
    }
  }

  // Compare the list of ids in this Spec to all previous Specs to ensure
  // uniqueness.
  Teuchos::RCP<UniqueSpecList> other_specs = unique_specs_[kind];
  if (other_specs == Teuchos::null) {
    other_specs = Teuchos::rcp(new UniqueSpecList());
    unique_specs_[kind] = other_specs;
  } else {
    for (auto uspec = other_specs->begin(); uspec != other_specs->end();
         ++uspec) {
      MeshIDs overlap;
      MeshIDs::iterator overlap_end;
      const MeshIDs& prev_spec_ids = *(*uspec)->second;

      std::set_intersection(prev_spec_ids.begin(),
                            prev_spec_ids.end(),
                            this_spec_ids->begin(),
                            this_spec_ids->end(),
                            std::inserter(overlap, overlap.end()));
      if (overlap.size() != 0) {
        Errors::Message m;
        m << "Conflicting definitions (overlapping) of id sets.";
        Exceptions::amanzi_throw(m);
      }
    }
  }

  // If we've gotten this far, things are OK. Add the spec to the lists.
  spec_list_.push_back(spec);
  other_specs->push_back(Teuchos::rcp(new UniqueSpec(spec, this_spec_ids)));
};

} // namespace Functions
} // namespace Amanzi
