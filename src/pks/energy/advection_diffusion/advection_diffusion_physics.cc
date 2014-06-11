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
  Teuchos::RCP<const CompositeVector> temp0 =
    S_inter_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> temp1 =
    S_next_->GetFieldData("temperature");

  Teuchos::RCP<const CompositeVector> cv0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cv1 =
    S_next_->GetFieldData("cell_volume");

  double dt = S_next_->time() - S_inter_->time();

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
  const Epetra_MultiVector& flux_f = *darcy_flux->ViewComponent("face",true);
  advection_->set_flux(darcy_flux);

  // put the advected quantity in cells
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  const Epetra_MultiVector& temp_f = *temp->ViewComponent("face", true);

  *field->ViewComponent("cell", false) = *temp->ViewComponent("cell", false);
  Epetra_MultiVector& field_f = *field->ViewComponent("face",true);

  // put the boundary fluxes in faces -- assumes all Dirichlet BC in temperature!
  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    field_f[0][f] = temp_f[0][f] * std::abs(flux_f[0][f]);
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
  // compute the stiffness matrix at the new time
  Teuchos::RCP<const CompositeVector> temp =
    S->GetFieldData("temperature");

  // get conductivity, and push it into whetstone tensor
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity");

  // calculate the div-grad operator, apply it to temperature, and add to residual
  matrix_->CreateMFDstiffnessMatrices(thermal_conductivity.ptr());
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->ComputeNegativeResidual(*temp, g.ptr());
};

} //namespace Energy
} //namespace Amanzi
