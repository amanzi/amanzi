/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "overland.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void OverlandFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(S_next_.ptr());

  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity", name_);
  matrix_->CreateMFDstiffnessMatrices(cond.ptr());
  matrix_->CreateMFDrhsVectors();

  // update the potential
  S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  if (update_flux_ == UPDATE_FLUX_ITERATION) {
    Teuchos::RCP<CompositeVector> flux =
        S->GetFieldData("surface_flux", name_);
    matrix_->DeriveFlux(*pres_elev, flux.ptr());
  }

  // Patch up BCs for zero-gradient
  FixBCsForOperator_(S_next_.ptr());

  // assemble the stiffness matrix
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, g.ptr());
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void OverlandFlow::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // get these fields
  Teuchos::RCP<const CompositeVector> h1 = S_next_->GetFieldData(key_);
  Teuchos::RCP<const CompositeVector> h0 = S_inter_->GetFieldData(key_);
  Teuchos::RCP<const CompositeVector> cv1 =
    S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP<const CompositeVector> cv0 =
    S_inter_->GetFieldData("surface_cell_volume");

  // Water content only has cells, while the residual has cells and faces.
  //  --   g <-- g - (cv*h)_t0/dt
  g->ViewComponent("cell",false)->Multiply(-1./dt,
          *cv0->ViewComponent("cell",false), *h0->ViewComponent("cell",false), 1.);
  //  --   g <-- g + (cv*h)_t1/dt
  g->ViewComponent("cell",false)->Multiply(1./dt,
          *cv1->ViewComponent("cell",false), *h1->ViewComponent("cell",false), 1.);

};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g) {
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  const Epetra_MultiVector& cv1 =
      *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);
  const Epetra_MultiVector& cv0 =
      *S_inter_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);

  if (is_source_term_) {
    // Add in external source term.
    S_next_->GetFieldEvaluator("surface_mass_source")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("surface_mass_source")
        ->HasFieldChanged(S_inter_.ptr(), name_);

    const Epetra_MultiVector& source0 =
        *S_inter_->GetFieldData("surface_mass_source")->ViewComponent("cell",false);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData("surface_mass_source")->ViewComponent("cell",false);

    //  --   g <-- g - 0.5 * (cv*h)_t0
    g->ViewComponent("cell",false)->Multiply(-0.5, cv0, source0, 1.);
    //  --   g <-- g - 0.5 * (cv*h)_t1
    g->ViewComponent("cell",false)->Multiply(-0.5, cv1, source1, 1.);
  }
};


} //namespace
} //namespace
