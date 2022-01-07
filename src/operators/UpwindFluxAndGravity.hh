/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Upwind a cell-centered field (e.g. rel perm) using first a given 
  face-based flux (e.g. Darcy flux), and then gravity.
*/

#ifndef AMANZI_UPWIND_FLUX_AND_GRAVITY_HH_
#define AMANZI_UPWIND_FLUX_AND_GRAVITY_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Mesh.hh"

// Operators
#include "Upwind.hh"
#include "UpwindFlux.hh"
#include "UpwindGravity.hh"

namespace Amanzi {
namespace Operators {

template<class Model>
class UpwindFluxAndGravity : public Upwind<Model> {
 public:
  UpwindFluxAndGravity(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                       Teuchos::RCP<const Model> model) :
      Upwind<Model>(mesh, model),
      upwind_flux_(mesh, model),
      upwind_gravity_(mesh, model) {};
  ~UpwindFluxAndGravity() {};

  // main methods
  // -- initialization of control parameters
  void Init(Teuchos::ParameterList& plist);

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly. 
  void Compute(const CompositeVector& flux, const CompositeVector& solution,
               const std::vector<int>& bc_model, CompositeVector& field);

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly. 
  Teuchos::RCP<CompositeVectorSpace> Map() {
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true)
       ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
       ->AddComponent("dirichlet_faces", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1)
       ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1)
       ->AddComponent("grav", AmanziMesh::Entity_kind::FACE, 1);
    return cvs;
  }

 private:
  using Upwind<Model>::mesh_;
  using Upwind<Model>::model_;
  using Upwind<Model>::face_comp_;

 private:
  int method_;
  UpwindFlux<Model> upwind_flux_;
  UpwindGravity<Model> upwind_gravity_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
template<class Model>
void UpwindFluxAndGravity<Model>::Init(Teuchos::ParameterList& plist)
{
  upwind_flux_.Init(plist);
  upwind_gravity_.Init(plist);

  method_ = Operators::OPERATOR_UPWIND_FLUX + Operators::OPERATOR_UPWIND_GRAVITY;
}


/* ******************************************************************
* Upwind field is placed in component "face" of field.
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
template<class Model>
void UpwindFluxAndGravity<Model>::Compute(
    const CompositeVector& flux, const CompositeVector& solution,
    const std::vector<int>& bc_model, CompositeVector& field)
{
  upwind_flux_.set_face_comp("face");
  upwind_flux_.Compute(flux, solution, bc_model, field);

  upwind_gravity_.set_face_comp("grav");
  upwind_gravity_.Compute(flux, solution, bc_model, field);
}

}  // namespace Operators
}  // namespace Amanzi

#endif

