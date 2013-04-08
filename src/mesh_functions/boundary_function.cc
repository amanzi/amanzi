/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon

Function applied to a mesh component with at most one function application per
entity.

------------------------------------------------------------------------- */

#include "boundary_function.hh"

namespace Amanzi {
namespace Functions {

void BoundaryFunction::Define(const std::vector<std::string> &regions,
        const Teuchos::RCP<const Function> &f) {

  // Create the domain
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));

  // add the spec
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
};


void BoundaryFunction::Define(std::string region,
        const Teuchos::RCP<const Function> &f) {

  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
};

// Evaluate values at time.
void BoundaryFunction::Compute(double time) {
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  if (specs_and_ids_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  double *args = new double[1+dim];
  double *xargs = args+1;
  args[0] = time;

  // Loop over all FACE specs and evaluate the function at all IDs in the
  // list.
  for (SpecAndIDsList::const_iterator
         spec_and_ids=specs_and_ids_[AmanziMesh::FACE]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::FACE]->end(); ++spec_and_ids) {
    // Here we could specialize on the argument signature of the function:
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. Right now we just assume
    // the most general case.
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->face_centroid(*id);
      for (int i=0; i!=dim; ++i) xargs[i] = xc[i];
      // Careful tracing of the typedefs is required here: spec_and_ids->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
    }
  }

  delete [] args;
};


void BoundaryFunction::Finalize() {
  finalized_ = true;
  if (specs_and_ids_.size() == 0) { return; }

  // Create the map of values, for now just setting up memory.
  for (SpecAndIDsList::const_iterator spec_and_ids =
         specs_and_ids_[AmanziMesh::FACE]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::FACE]->end();
       ++spec_and_ids) {
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id=ids->begin(); id!=ids->end(); ++id) {
      value_[*id];
    };
  }

  //TODO: Verify that the faces in this_domain are all boundary faces.
};

} // namespace
} // namespace
