/*
  Transport PK

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#include "TransportBoundaryFunction_Tracer.hh"

namespace Amanzi {
namespace Transport {

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
  if (unique_specs_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = time;

  // Loop over side set specs and evaluate the function at all faces 
  // in the side set list.
  int n = 0;
  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::FACE]->begin();
       uspec != unique_specs_[AmanziMesh::FACE]->end(); ++uspec) {
    // We could specialize on the argument signature of the function:
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. 
    // Right now we just assume the most general case.
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->face_centroid(*id);
      for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
      values_[n++][0] = (*(*uspec)->first->second)(args)[0];
    }
  }
}


/* ******************************************************************
* Generate space for data (face ids and values).
****************************************************************** */
void TransportBoundaryFunction_Tracer::Finalize_() 
{
  if (unique_specs_.size() == 0) return;

  std::vector<double> v;
  v.push_back(0.0);

  int n = 0;
  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::FACE]->begin();
       uspec != unique_specs_[AmanziMesh::FACE]->end(); ++uspec) {
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      faces_.push_back(*id);
      values_.push_back(v);
      n++;
    }
  }
}

}  // namespace Transport
}  // namespace Amanzi

