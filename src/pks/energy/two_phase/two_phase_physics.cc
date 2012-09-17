/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Solves:

du/dt + v dot grad h = div Ke grad T .......... fix me!
------------------------------------------------------------------------- */

#include "advection.hh"
#include "eos.hh"
#include "iem.hh"
#include "field_evaluator.hh"
#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Accumulation of internal energy term du/dt
// -------------------------------------------------------------
void TwoPhase::AddAccumulation_(Teuchos::RCP<CompositeVector> g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the energy at both the old and new times.
  S_next_->GetFieldEvaluator("energy")->HasFieldChanged(S_next_.ptr(), "energy_pk");
  S_inter_->GetFieldEvaluator("energy")->HasFieldChanged(S_inter_.ptr(), "energy_pk");

  // get the energy at each time
  Teuchos::RCP<const CompositeVector> e1 = S_next_->GetFieldData("energy");
  Teuchos::RCP<const CompositeVector> e0 = S_inter_->GetFieldData("energy");

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  g->ViewComponent("cell", false)
    ->Update(1.0/dt, *e1->ViewComponent("cell", false),
          -1.0/dt, *e0->ViewComponent("cell", false), 1.0);
};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, v dot grad h.
// -------------------------------------------------------------
void TwoPhase::AddAdvection_(const Teuchos::RCP<State> S,
        const Teuchos::RCP<CompositeVector> g, bool negate) {

  Teuchos::RCP<CompositeVector> field = advection_->field();
  field->PutScalar(0);

  // set the flux field as the darcy flux
  // NOTE: darcy_flux is a MOLAR flux by choice of the richards flow pk, i.e.
  // [flux] =  mol/(m^2*s)

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no FieldEvaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Richards::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetFieldEvaluator("darcy_flux")->HasFieldChanged(S.ptr(), "energy_pk");
  Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux");
  advection_->set_flux(darcy_flux);

  // put the advected quantity in cells
  S->GetFieldEvaluator("enthalpy_liquid")->HasFieldChanged(S.ptr(), "energy_pk");
  Teuchos::RCP<const CompositeVector> enthalpy = S->GetFieldData("enthalpy_liquid");
  *field->ViewComponent("cell", false) = *enthalpy->ViewComponent("cell", false);


  // put the boundary fluxes in faces -- assumes all Dirichlet BC in temperature!
  // NOTE this boundary flux is in enthalpy, and
  // h = n(T,p) * u_l(T) + p_l
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");

  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    double p = (*pres)("face",f);
    double T = bc->second;
    double dens = eos_liquid_->MolarDensity(T,p);
    double int_energy = iem_liquid_->InternalEnergy(T);
    double enthalpy = int_energy + p/dens;

    (*field)("face",f) = enthalpy * fabs((*darcy_flux)("face",f));
  }

  // apply the advection operator and add to residual
  advection_->Apply(bc_flux_);
  if (negate) {
    for (int c=0; c!=g->size("cell"); ++c) {
      (*g)("cell",c) -= (*field)("cell",c);
    }
  } else {
    for (int c=0; c!=g->size("cell"); ++c) {
      (*g)("cell",c) += (*field)("cell",c);
    }
  }
};


// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void TwoPhase::ApplyDiffusion_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> g) {
  // update the thermal conductivity
  S->GetFieldEvaluator("thermal_conductivity")->HasFieldChanged(S.ptr(), "energy_pk");
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(*thermal_conductivity);
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  matrix_->ComputeNegativeResidual(*temp, g);
};

} //namespace Energy
} //namespace Amanzi
