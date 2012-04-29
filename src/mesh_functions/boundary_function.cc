#include "boundary_function.hh"

#include <algorithm>
#include "errors.hh"

namespace Amanzi {

void BoundaryFunction::Define(const std::vector<std::string> &regions,
                              const Teuchos::RCP<const Function> &f)
{
  // Form the set of face indices that belong to any of the given regions.
  // We use a std::set here to filter out any duplicate indices.
  std::set<AmanziMesh::Entity_ID> this_domain;
  for (std::vector<std::string>::const_iterator r = regions.begin(); r != regions.end(); ++r) {
    if ((*mesh_).valid_set_name(*r, AmanziMesh::FACE)) {
      AmanziMesh::Entity_ID_List face_list;
      (*mesh_).get_set_entities(*r, AmanziMesh::FACE, AmanziMesh::USED, &face_list);
      this_domain.insert(face_list.begin(), face_list.end());
    } else {
      Errors::Message m;
      m << "unknown region: \"" << r->c_str() << "\"";
      Exceptions::amanzi_throw(m);
    }
  }
  
  //TODO: Verify that the faces in this_domain are all boundary faces.
  
  // Verify that the faces in this_domain are disjoint from the previous domains.
  for (SpecList::const_iterator s = spec_list_.begin(); s != spec_list_.end(); ++s) {
    std::set<AmanziMesh::Entity_ID> overlap;
    std::set<AmanziMesh::Entity_ID>::iterator overlap_end;
    const std::set<AmanziMesh::Entity_ID> &prev_domain = s->first;
    std::set_intersection(prev_domain.begin(), prev_domain.end(),
                          this_domain.begin(), this_domain.end(),
                          std::inserter(overlap, overlap.end()));
    if (overlap.size() != 0) {
      Errors::Message m;
      m << "conflicting definition";
      Exceptions::amanzi_throw(m);
    }
  }
  
  // Add this specification to the list.
  spec_list_.push_back(std::make_pair(this_domain, f));
  
  // Go ahead and register the faces in the value_ map; this enables one to
  // get at the faces in the boundary function without computing its value.
  for (Domain::const_iterator d = this_domain.begin(); d != this_domain.end(); ++d) value_[*d];
}

} // namespace Amanzi
