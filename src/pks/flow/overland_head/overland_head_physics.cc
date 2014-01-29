/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "overland_head.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void OverlandHeadFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
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
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, g.ptr());


  std::cout << std::setprecision(15) << "FUCK RESIDUAL F(237): h+z, r = " << (*pres_elev->ViewComponent("face",false))[0][237] << ", " 
            << (*g->ViewComponent("face",false))[0][237] << std::endl;
  std::cout << std::setprecision(15) << "            C(78,79): h+z = " << (*pres_elev->ViewComponent("cell",false))[0][78] << ", "
            << (*pres_elev->ViewComponent("cell",false))[0][79] << std::endl;
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void OverlandHeadFlow::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // get these fields
  S_next_->GetFieldEvaluator("surface_water_content")
      ->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator("surface_water_content")
      ->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> wc1 =
      S_next_->GetFieldData("surface_water_content");
  Teuchos::RCP<const CompositeVector> wc0 =
      S_inter_->GetFieldData("surface_water_content");

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell",false)->Update(1.0/dt, *wc1->ViewComponent("cell",false),
          -1.0/dt, *wc0->ViewComponent("cell",false), 1.0);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandHeadFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g) {
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  const Epetra_MultiVector& cv1 =
      *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);

  if (is_source_term_) {
    // Add in external source term.
    S_next_->GetFieldEvaluator("surface_mass_source")
        ->HasFieldChanged(S_next_.ptr(), name_);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData("surface_mass_source")->ViewComponent("cell",false);

    // External source term is in [m water / s], not in [mols / s], so a
    // density is required.  This density should be upwinded.
    S_next_->GetFieldEvaluator("surface_molar_density_liquid")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_next_->GetFieldEvaluator("surface_source_molar_density")
        ->HasFieldChanged(S_next_.ptr(), name_);

    const Epetra_MultiVector& nliq1 =
        *S_next_->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& nliq1_s =
        *S_next_->GetFieldData("surface_source_molar_density")
        ->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      double s1 = source1[0][c] > 0. ? source1[0][c] * nliq1_s[0][c] : source1[0][c] * nliq1[0][c];
      g_c[0][c] -= cv1[0][c] * s1;
    }
  }

  if (coupled_to_subsurface_via_head_) {
    // Add in source term from coupling.
    S_next_->GetFieldEvaluator("overland_source_from_subsurface")
        ->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> source1 =
        S_next_->GetFieldData("overland_source_from_subsurface");

    // source term is in units of [mol / s]
    g_c.Update(-1., *source1->ViewComponent("cell",false), 1.);
  }
};


} //namespace
} //namespace
