/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

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

class UpwindFluxAndGravity : public Upwind {
 public:
  UpwindFluxAndGravity(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Upwind(mesh), upwind_flux_(mesh), upwind_gravity_(mesh){};
  ~UpwindFluxAndGravity(){};

  // main methods
  // -- initialization of control parameters
  void Init(Teuchos::ParameterList& plist);

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly.
  void
  Compute(const CompositeVector& flux, const std::vector<int>& bc_model, CompositeVector& field);

  // -- returns combined map for the original and upwinded fields.
  // -- Currently, composite vector cannot be extended on a fly.
  Teuchos::RCP<CompositeVectorSpace> Map()
  {
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1)
      ->AddComponent("grav", AmanziMesh::FACE, 1);
    return cvs;
  }

 private:
  int method_;
  UpwindFlux upwind_flux_;
  UpwindGravity upwind_gravity_;
};


/* ******************************************************************
* Public init method. It is not yet used.
****************************************************************** */
inline void
UpwindFluxAndGravity::Init(Teuchos::ParameterList& plist)
{
  upwind_flux_.Init(plist);
  upwind_gravity_.Init(plist);

  method_ = Operators::OPERATOR_UPWIND_FLUX + Operators::OPERATOR_UPWIND_GRAVITY;
}


/* ******************************************************************
* Upwind field is placed in component "face" of field.
* Upwinded field must be calculated on all faces of the owned cells.
****************************************************************** */
inline void
UpwindFluxAndGravity::Compute(const CompositeVector& flux,
                              const std::vector<int>& bc_model,
                              CompositeVector& field)
{
  upwind_flux_.set_face_comp("face");
  upwind_flux_.Compute(flux, bc_model, field);

  upwind_gravity_.set_face_comp("grav");
  upwind_gravity_.Compute(flux, bc_model, field);
}

} // namespace Operators
} // namespace Amanzi

#endif
