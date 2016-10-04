/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void OverlandPressureFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(S_next_.ptr());

  // update the stiffness matrix
  matrix_->Init();
  Teuchos::RCP<const CompositeVector> cond =
<<<<<<< HEAD
    S_next_->GetFieldData(getKey(domain_,"upwind_overland_conductivity"), name_);
=======
    S_next_->GetFieldData("upwind_overland_conductivity", name_);
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update the potential
<<<<<<< HEAD
  S->GetFieldEvaluator(getKey(domain_,"pres_elev"))->HasFieldChanged(S.ptr(), name_);
=======
  S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  // Patch up BCs for zero-gradient
  FixBCsForOperator_(S_next_.ptr());

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
<<<<<<< HEAD
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData(getKey(domain_,"pres_elev"));
  if (update_flux_ == UPDATE_FLUX_ITERATION) {
    Teuchos::RCP<CompositeVector> flux =
      S->GetFieldData(getKey(domain_,"mass_flux"), name_);
=======
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  if (update_flux_ == UPDATE_FLUX_ITERATION) {
    Teuchos::RCP<CompositeVector> flux =
        S->GetFieldData("surface-mass_flux", name_);
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
    matrix_diff_->UpdateFlux(*pres_elev, *flux);
  }

  // assemble the stiffness matrix
  matrix_diff_->ApplyBCs(true, true);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, *g);
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void OverlandPressureFlow::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // get these fields
<<<<<<< HEAD
  S_next_->GetFieldEvaluator(conserved_key_)
      ->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(conserved_key_)
      ->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> wc1 =
      S_next_->GetFieldData(conserved_key_);
  Teuchos::RCP<const CompositeVector> wc0 =
      S_inter_->GetFieldData(conserved_key_);
=======
  S_next_->GetFieldEvaluator("surface-water_content")
      ->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator("surface-water_content")
      ->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> wc1 =
      S_next_->GetFieldData("surface-water_content");
  Teuchos::RCP<const CompositeVector> wc0 =
      S_inter_->GetFieldData("surface-water_content");
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell",false)->Update(1.0/dt, *wc1->ViewComponent("cell",false),
          -1.0/dt, *wc0->ViewComponent("cell",false), 1.0);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandPressureFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g) {
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  const Epetra_MultiVector& cv1 =
<<<<<<< HEAD
    *S_next_->GetFieldData(getKey(domain_,"cell_volume"))->ViewComponent("cell",false);
=======
      *S_next_->GetFieldData("surface-cell_volume")->ViewComponent("cell",false);
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  if (is_source_term_) {
    // Add in external source term.
    S_next_->GetFieldEvaluator(mass_source_key_)
        ->HasFieldChanged(S_next_.ptr(), name_);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData(mass_source_key_)->ViewComponent("cell",false);

    if (source_in_meters_) {
      // External source term is in [m water / s], not in [mols / s], so a
      // density is required.  This density should be upwinded.
<<<<<<< HEAD
      S_next_->GetFieldEvaluator(getKey(domain_,"molar_density_liquid"))
          ->HasFieldChanged(S_next_.ptr(), name_);
      S_next_->GetFieldEvaluator(getKey(domain_,"source_molar_density"))
          ->HasFieldChanged(S_next_.ptr(), name_);

      const Epetra_MultiVector& nliq1 =
        *S_next_->GetFieldData(getKey(domain_,"molar_density_liquid"))
          ->ViewComponent("cell",false);
      const Epetra_MultiVector& nliq1_s =
        *S_next_->GetFieldData(getKey(domain_,"source_molar_density"))
=======
      S_next_->GetFieldEvaluator("surface-molar_density_liquid")
          ->HasFieldChanged(S_next_.ptr(), name_);
      S_next_->GetFieldEvaluator("surface-source_molar_density")
          ->HasFieldChanged(S_next_.ptr(), name_);

      const Epetra_MultiVector& nliq1 =
          *S_next_->GetFieldData("surface-molar_density_liquid")
          ->ViewComponent("cell",false);
      const Epetra_MultiVector& nliq1_s =
        *S_next_->GetFieldData("surface-source_molar_density")
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
          ->ViewComponent("cell",false);

      int ncells = g_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        double s1 = source1[0][c] > 0. ? source1[0][c] * nliq1_s[0][c] : source1[0][c] * nliq1[0][c];
        g_c[0][c] -= cv1[0][c] * s1;
      }
    } else {
      int ncells = g_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        g_c[0][c] -= cv1[0][c] * source1[0][c];
      }
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
<<<<<<< HEAD
  /*
  //I-ADDED
  if (coupled_to_subsurface_via_flux_) {
    Epetra_MultiVector& source = *S_next_->GetFieldData(getKey(domain_.substr(0,domain_.size()-8),"surface_subsurface_flux"),
                                                        S_inter_->GetField(getKey(domain_.substr(0,domain_.size()-8),"surface_subsurface_flux"))->owner())->ViewComponent("cell",false);
    source = g_c;
    g_c =  *g->ViewComponent("cell",false)->PutScalar(0.0);
  }
  */
=======
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
};


} //namespace
} //namespace
