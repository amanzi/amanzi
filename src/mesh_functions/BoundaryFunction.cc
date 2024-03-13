/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: (v1) Neil Carlson
           (v2) Ethan Coon
*/

/*
  Mesh Functions

  Function applied to a mesh component with at most one function
  application per entity.
*/

#include "BoundaryFunction.hh"

namespace Amanzi {
namespace Functions {

void
BoundaryFunction::Define(const std::vector<std::string>& regions,
                         const Teuchos::RCP<const MultiFunction>& f)
{
  // Create the domain
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::Entity_kind::FACE));

  // add the spec
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
};


void
BoundaryFunction::Define(std::string& region, const Teuchos::RCP<const MultiFunction>& f)
{
  RegionList regions(1, region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::Entity_kind::FACE));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
};


// Evaluate values at time.
void
BoundaryFunction::Compute(double time)
{
  // lazily generate space for the values
  if (!finalized_) { Finalize(); }

  if (unique_specs_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->getSpaceDimension();
  std::vector<double> args(1 + dim);
  args[0] = time;

  // Loop over all FACE specs and evaluate the function at all IDs in the
  // list.
  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::Entity_kind::FACE]->begin();
       uspec != unique_specs_[AmanziMesh::Entity_kind::FACE]->end();
       ++uspec) {
    // Here we could specialize on the argument signature of the function:
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. Right now we just assume
    // the most general case.
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->getFaceCentroid(*id);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // Careful tracing of the typedefs is required here: uspec->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*uspec)->first->second)(args)[0];
    }
  }
};


void
BoundaryFunction::Finalize()
{
  finalized_ = true;
  if (unique_specs_.size() == 0) { return; }

  // Create the map of values, for now just setting up memory.
  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::Entity_kind::FACE]->begin();
       uspec != unique_specs_[AmanziMesh::Entity_kind::FACE]->end();
       ++uspec) {
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) { value_[*id]; };
  }

  //TODO: Verify that the faces in this_domain are all boundary faces.
};

} // namespace Functions
} // namespace Amanzi
