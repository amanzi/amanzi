#include "boundary_function.hh"

#include <algorithm>
#include "errors.hh"

namespace Amanzi {

/* ****************************************************************
* Populates internal array with function values.
**************************************************************** */
void BoundaryFunction::Compute(double t)
{
  int dim = (*mesh_).space_dimension();
  double *args = new double[1+dim];
  double *xargs = args+1;
  args[0] = t;

  for (SpecList::const_iterator s = spec_list_.begin(); s != spec_list_.end(); ++s) {
    // Here we could specialize on the argument signature of the function: 
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. Right now we just assume
    // the most general case.
    const Domain& domain = s->first;
    for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
      const AmanziGeometry::Point& xc = mesh_->face_centroid(*d);
      for (int i = 0; i < dim; ++i) xargs[i] = xc[i];
      value_[*d] = (*(s->second))(args);
    }
  }
  delete [] args;
}


/* ****************************************************************
* Populates internal array with function values. 
**************************************************************** */
void BoundaryFunction::Define(const std::vector<std::string> &regions,
                              const Teuchos::RCP<const Function> &f)
{
  // Form the set of face indices that belong to any of the given regions.
  // We use a std::set here to filter out any duplicate indices.
  std::set<AmanziMesh::Entity_ID> this_domain;
  for (std::vector<std::string>::const_iterator r = regions.begin(); r != regions.end(); ++r) {
    if ((*mesh_).valid_set_name(*r, AmanziMesh::FACE)) {
      AmanziMesh::Entity_ID_List face_list;
      (*mesh_).get_set_entities(*r, AmanziMesh::FACE, AmanziMesh::OWNED, &face_list);
      this_domain.insert(face_list.begin(), face_list.end());
    } else {
      Errors::Message m;
      m << "\n   unknown region name or wrong topological dimension for \"" << r->c_str() << "\"";
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
