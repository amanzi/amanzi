/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Mesh.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void OverlandFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
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
        S->GetFieldData("surface_flux", name_);
    matrix_->DeriveFlux(*pres_elev, flux.ptr());
    flux->ScatterMasterToGhosted();
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
  Teuchos::RCP <const CompositeVector> cv1 =
      S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP <const CompositeVector> cv0 =
      S_inter_->GetFieldData("surface_cell_volume");

  if (is_source_term_) {
    S_next_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_inter_.ptr(), name_);

    Teuchos::RCP <const CompositeVector> source0 =
        S_inter_->GetFieldData("overland_source");
    Teuchos::RCP <const CompositeVector> source1 =
        S_next_->GetFieldData("overland_source");

    //  --   g <-- g - 0.5 * (cv*h)_t0
    g->ViewComponent("cell",false)->Multiply(-0.5,
            *cv0->ViewComponent("cell",false),
            *source0->ViewComponent("cell",false), 1.);
    //  --   g <-- g - 0.5 * (cv*h)_t1
    g->ViewComponent("cell",false)->Multiply(-0.5,
            *cv1->ViewComponent("cell",false),
            *source1->ViewComponent("cell",false), 1.);
  }

  if (is_coupling_term_) {
    S_next_->GetFieldEvaluator("overland_source_from_subsurface")
        ->HasFieldChanged(S_next_.ptr(), name_);

    Teuchos::RCP<const CompositeVector> pres = S_next_ ->GetFieldData(key_);
    Teuchos::RCP<const CompositeVector> source1 =
        S_next_->GetFieldData("overland_source_from_subsurface");

    //  --   g <-- g - (cv*h)_t1
    g->ViewComponent("cell",false)->Multiply(-1.,
            *cv1->ViewComponent("cell",false),
            *source1->ViewComponent("cell",false), 1.);
  }
};


} //namespace
} //namespace
