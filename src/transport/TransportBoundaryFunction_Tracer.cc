/*
  This is the transport component of the Amanzi code. 

  License: see $AMANZI_DIR/COPYRIGHT
  Author (v1): Neil Carlson
         (v2): Ethan Coon
         (v3): Konstantin Lipnikov
*/

#include "TransportBoundaryFunction_Tracer.hh"

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Tracer::Define(
    const std::vector<std::string> &regions,
    const Teuchos::RCP<const MultiFunction> &f) 
{
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  Finalize_();
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Tracer::Define(
    std::string region, const Teuchos::RCP<const MultiFunction> &f) 
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  Finalize_();
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportBoundaryFunction_Tracer::Compute(double time) {
  if (specs_and_ids_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  double *args = new double[dim + 1];
  double *xargs = args + 1;
  args[0] = time;

  // Loop over side set specs and evaluate the function at all faces 
  // in the side set list.
  int n = 0;
  for (SpecAndIDsList::const_iterator
       spec_and_ids = specs_and_ids_[AmanziMesh::FACE]->begin();
       spec_and_ids != specs_and_ids_[AmanziMesh::FACE]->end(); ++spec_and_ids) {
    // We could specialize on the argument signature of the function:
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. 
    // Right now we just assume the most general case.
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->face_centroid(*id);
      for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
      values_[n++][0] = (*(*spec_and_ids)->first->second)(args)[0];
    }
  }
  delete [] args;
}


/* ******************************************************************
* Generate space for data (face ids and values).
****************************************************************** */
void TransportBoundaryFunction_Tracer::Finalize_() 
{
  if (specs_and_ids_.size() == 0) return;

  std::vector<double> v;
  v.push_back(0.0);

  int n = 0;
  for (SpecAndIDsList::const_iterator spec_and_ids =
         specs_and_ids_[AmanziMesh::FACE]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::FACE]->end();
       ++spec_and_ids) {
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id=ids->begin(); id!=ids->end(); ++id) {
      faces_.push_back(*id);
      values_.push_back(v);
      n++;
    }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi

