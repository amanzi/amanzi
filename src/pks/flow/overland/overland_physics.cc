/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Mesh.hh"
#include "Mesh_MSTK.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void OverlandFlow::ApplyDiffusion_(const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<CompositeVector>& g) {
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
    matrix_->DeriveFlux(*pres_elev, flux);
  }

  // assemble the stiffness matrix
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, g);
};


// -------------------------------------------------------------
// Accumulation of internal energy term du/dt
// -------------------------------------------------------------
void OverlandFlow::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
  const CompositeVector & pres0        = *(S_inter_->GetFieldData(key_));
  const CompositeVector & pres1        = *(S_next_ ->GetFieldData(key_));
  const CompositeVector & cell_volume0 = *(S_inter_->GetFieldData("surface_cell_volume"));
  const CompositeVector & cell_volume1 = *(S_next_ ->GetFieldData("surface_cell_volume"));

  double dt = S_next_->time() - S_inter_->time();

  int c_owned = g->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",0,c) += (cell_volume1("cell",c)*pres1("cell",c)
                         - cell_volume0("cell",c)*pres0("cell",c))/dt ;
  }
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandFlow::AddLoadValue_(const Teuchos::RCP<CompositeVector>& g) {
  Teuchos::RCP <const CompositeVector> cell_volume1 =
      S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP <const CompositeVector> cell_volume0 =
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

    int c_owned = g->size("cell");
    for (int c=0; c!=c_owned; ++c) {
      (*g)("cell",c) -= ( (*source1)("cell",c) * (*cell_volume1)("cell",c)
                          + (*source0)("cell",c) * (*cell_volume0)("cell",c)) / 2.0;
    }
  }

  if (is_coupling_term_) {
    S_next_->GetFieldEvaluator("overland_source_from_subsurface")
        ->HasFieldChanged(S_next_.ptr(), name_);

    Teuchos::RCP<const CompositeVector> pres = S_next_ ->GetFieldData(key_);
    Teuchos::RCP<const CompositeVector> source1 =
        S_next_->GetFieldData("overland_source_from_subsurface");

    int c_owned = g->size("cell");
    for (int c=0; c!=c_owned; ++c) {
      //      (*g)("cell",c) -= std::max<double>((*source1)("cell",c), -abs((*pres)("cell",c)) / (S_next_->time() - S_inter_->time())) * (*cell_volume1)("cell",c);
      (*g)("cell",c) -= (*source1)("cell",c) * (*cell_volume1)("cell",c);
    }
  }
};


} //namespace
} //namespace
