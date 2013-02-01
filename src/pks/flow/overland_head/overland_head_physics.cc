/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Mesh.hh"
#include "Mesh_MSTK.hh"

#include "overland_head.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void OverlandHeadFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  // update the rel perm according to the scheme of choice
  bool update = UpdatePermeabilityData_(S.ptr());

  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> cond =
    S->GetFieldData("upwind_overland_conductivity", name_);
  matrix_->CreateMFDstiffnessMatrices(cond.ptr());

  // update the potential
  update |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  if (update_flux_ == UPDATE_FLUX_ITERATION) {
    Teuchos::RCP<CompositeVector> flux =
        S->GetFieldData("overland_flux", name_);
    pres_elev->ScatterMasterToGhosted();
    matrix_->DeriveFlux(*pres_elev, flux.ptr());
  }

  // assemble the stiffness matrix
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, g.ptr());
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
  Teuchos::RCP<const CompositeVector> cv1 =
    S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP<const CompositeVector> cv0 =
    S_inter_->GetFieldData("surface_cell_volume");

  // Water content only has cells, while the residual has cells and faces.
  //  --   g <-- g - (cv*h)_t0/dt
  g->ViewComponent("cell",false)->Multiply(-1./dt,
          *cv0->ViewComponent("cell",false), *wc0->ViewComponent("cell",false), 1.);
  //  --   g <-- g + (cv*h)_t1/dt
  g->ViewComponent("cell",false)->Multiply(1./dt,
          *cv1->ViewComponent("cell",false), *wc1->ViewComponent("cell",false), 1.);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandHeadFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g) {
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  const Epetra_MultiVector& cv1 =
      *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);
  const Epetra_MultiVector& cv0 =
      *S_inter_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);

  if (is_source_term_) {
    S_next_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_inter_.ptr(), name_);

    const Epetra_MultiVector& source0 =
        *S_inter_->GetFieldData("overland_source")->ViewComponent("cell",false);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData("overland_source")->ViewComponent("cell",false);

    // source term is in m water / s, not in mols / s.
    S_next_->GetFieldEvaluator("surface_molar_density_liquid")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("surface_molar_density_liquid")
        ->HasFieldChanged(S_inter_.ptr(), name_);

    const Epetra_MultiVector& nliq0 =
        *S_inter_->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& nliq1 =
        *S_next_->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);

    int ncells = g_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      g_c[0][c] -= 0.5* (cv0[0][c] * source0[0][c] * nliq0[0][c]
                         + cv1[0][c] * source1[0][c] * nliq1[0][c]);
    }
  }

  if (is_coupling_term_) {
    S_next_->GetFieldEvaluator("overland_source_from_subsurface")
        ->HasFieldChanged(S_next_.ptr(), name_);

    Teuchos::RCP<const CompositeVector> source1 =
        S_next_->GetFieldData("overland_source_from_subsurface");

    //  --   g <-- g - (cv*h)_t1
    // source term is in units of mol / m^2 / s, where the factor of face_area
    // / cell_volume has already been dealt with.
    //    g_c.Multiply(-1., cv1, *source1->ViewComponent("cell",false), 1.);

    // source term is in units of mol / s
    g_c.Update(-1., *source1->ViewComponent("cell",false), 1.);
  }
};


} //namespace
} //namespace
