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

  Teuchos::RCP<const CompositeVector> cell_volume0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cell_volume1 =
    S_next_->GetFieldData("cell_volume");

  double dt = S_next_->time() - S_inter_->time();

  // NOTE: gas and liquid are done in a ?? basis, but rock is done in a mass basis

  int c_owned = g->size("cell",false);
  for (int c=0; c != c_owned; ++c) {
    // add the time derivative of energy density to the residual
    (*g)("cell",c) += ((*temp1)("cell",c)*(*cell_volume1)("cell",c)
                         - (*temp0)("cell",c)*(*cell_volume0)("cell",c))/dt;
  }
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

  // put the boundary fluxes in faces -- assumes all Dirichlet BC in temperature!
  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    // note that temp has the correct BC value in its face already, as it was
    // solved from the previous timestep.
    (*field)("face",f) = (*temp)("face",f) * fabs((*darcy_flux)("face",f));
  }

  // apply the advection operator and add to residual
  advection_->Apply(bc_flux_);
  int c_owned = g->size("cell",false);
  if (negate) {
    for (int c=0; c!=c_owned; ++c) {
      (*g)("cell",c) -= (*field)("cell",c);
    }
  } else {
    for (int c=0; c!=c_owned; ++c) {
      (*g)("cell",c) += (*field)("cell",c);
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
  matrix_->AssembleGlobalMatrices();
  matrix_->ComputeNegativeResidual(*temp, g);
};

} //namespace Energy
} //namespace Amanzi
