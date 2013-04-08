/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Function applied to a mesh component with at most one function application per
entity.

------------------------------------------------------------------------- */

#include "errors.hh"

#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Functions {

// Overload the AddSpec method to check uniqueness.
void UniqueMeshFunction::AddSpec(const Teuchos::RCP<Spec>& spec) {
  // Ensure uniqueness of the spec and create the set of IDs contained in the
  // Domain of the spec.
  Teuchos::RCP<Domain> domain = spec->first;
  AmanziMesh::Entity_kind kind = domain->second;

  // Loop over regions in the spec, getting their ids and adding to the set.
  Teuchos::RCP<SpecIDs> this_spec_ids = Teuchos::rcp(new SpecIDs());
  for (RegionList::const_iterator region=domain->first.begin();
         region!=domain->first.end(); ++region) {

    // Get the ids from the mesh by region name and entity kind.
    if (mesh_->valid_set_name(*region, kind)) {
      AmanziMesh::Entity_ID_List id_list;
      mesh_->get_set_entities(*region, kind, AmanziMesh::USED, &id_list);
      this_spec_ids->insert(id_list.begin(), id_list.end());
    } else {
      std::stringstream m;
      m << "unknown region: \"" << *region << "\"";
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
    }
  }

  // Compare the list of ids in this Spec to all previous Specs to ensure
  // uniqueness.
  Teuchos::RCP<SpecAndIDsList> other_specs = specs_and_ids_[kind];
  if (other_specs == Teuchos::null) {
    other_specs = Teuchos::rcp(new SpecAndIDsList());
    specs_and_ids_[kind] = other_specs;
  } else {
    for (SpecAndIDsList::const_iterator spec_and_ids = other_specs->begin();
         spec_and_ids!=other_specs->end(); ++spec_and_ids) {

      std::set<AmanziMesh::Entity_ID> overlap;
      std::set<AmanziMesh::Entity_ID>::iterator overlap_end;
      const std::set<AmanziMesh::Entity_ID> &prev_spec_ids = *(*spec_and_ids)->second;
      std::set_intersection(prev_spec_ids.begin(), prev_spec_ids.end(),
                            this_spec_ids->begin(), this_spec_ids->end(),
                            std::inserter(overlap, overlap.end()));
      if (overlap.size() != 0) {
        Errors::Message m;
        m << "conflicting definition";
        Exceptions::amanzi_throw(m);
      }
    }
  }
  // If we've gotten this far, things are OK.  Add the spec to the lists.
  spec_list_.push_back(spec);
  other_specs->push_back(Teuchos::rcp(new SpecAndIDs(spec, this_spec_ids)));
};


} //namespace
} //namespace

