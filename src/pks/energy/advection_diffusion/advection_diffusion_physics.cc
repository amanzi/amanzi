/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Physics methods of the advection-diffusion equation.

Solves:

dT/dt  +  u dot grad(T) = div K grad T

------------------------------------------------------------------------- */

#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {


// dT/dt portion of the residual function
void AdvectionDiffusion::AddAccumulation_(Teuchos::RCP<CompositeVector> g) {
  S_next_->GetFieldEvaluator("temperature")->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator("temperature")->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> temp0 =
    S_inter_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> temp1 =
    S_next_->GetFieldData("temperature");

  Teuchos::RCP<const CompositeVector> cv0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cv1 =
    S_next_->GetFieldData("cell_volume");

  double dt = S_next_->time() - S_inter_->time();
  ASSERT(dt > 0.);

  //  --   g <-- g - (cv*h)_t0/dt
  g->ViewComponent("cell",false)->Multiply(-1./dt,
          *cv0->ViewComponent("cell",false), *temp0->ViewComponent("cell",false), 1.);
  //  --   g <-- g + (cv*h)_t1/dt
  g->ViewComponent("cell",false)->Multiply(1./dt,
          *cv1->ViewComponent("cell",false), *temp1->ViewComponent("cell",false), 1.);
};



// u dot grad T portion of the residual function
void AdvectionDiffusion::AddAdvection_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> g, bool negate) {

  Teuchos::RCP<CompositeVector> field = advection_->field();
  field->PutScalar(0);

  // set the flux field as the darcy flux
  Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux");
  advection_->set_flux(darcy_flux);

  // put the advected quantity in cells
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  *field->ViewComponent("cell", false) = *temp->ViewComponent("cell", false);


  // put the boundary fluxes in incoming faces for Dirichlet BCs.
  const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")
      ->ViewComponent("face",false);
  // put the boundary fluxes in faces -- assumes all Dirichlet BC in temperature!
  Epetra_MultiVector& field_f = *field->ViewComponent("face",true);
  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    //    std::cout << "putting Adv bc in face at " << mesh_->face_centroid(bc->first) << " at val " << bc->second << " with flux " << flux[0][bc->first] << std::endl;
    int f = bc->first;
    field_f[0][f] = bc->second * std::abs(flux[0][f]);
  }

  // apply the advection operator and add to residual
  advection_->Apply(bc_flux_);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);
  Teuchos::RCP<const CompositeVector> field_const(field);
  const Epetra_MultiVector& field_c =
      *field_const->ViewComponent("cell", false);

  int c_owned = g_c.MyLength();
  if (negate) {
    for (int c=0; c!=c_owned; ++c) {
      g_c[0][c] -= field_c[0][c];
    }
  } else {
    for (int c=0; c!=c_owned; ++c) {
      g_c[0][c] += field_c[0][c];
    }
  }
};


// - div K grad T part of the residual
void AdvectionDiffusion::ApplyDiffusion_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> g) {
  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity");

  matrix_->Init();
  matrix_diff_->Setup(thermal_conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // calculate the div-grad operator, apply it to temperature, and add to residual
  Teuchos::RCP<const CompositeVector> temp =
    S->GetFieldData("temperature");


  // finish assembly of the stiffness matrix
  matrix_diff_->ApplyBCs(true);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*temp, *g);
};

} //namespace Energy
} //namespace Amanzi
