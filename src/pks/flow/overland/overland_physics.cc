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

#if 0
}}
#endif

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void OverlandFlow::ApplyDiffusion_(const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<CompositeVector>& g) {
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  Teuchos::RCP<CompositeVector> darcy_flux =
    S->GetFieldData("overland_flux", "overland_flow");

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> cond =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(cond.ptr());
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
  const CompositeVector & pres0        = *(S_inter_->GetFieldData("overland_pressure"));
  const CompositeVector & pres1        = *(S_next_ ->GetFieldData("overland_pressure"));
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
        ->HasFieldChanged(S_next_.ptr(), "overland_pk");
    S_inter_->GetFieldEvaluator("overland_source")
        ->HasFieldChanged(S_inter_.ptr(), "overland_pk");

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
    S_next_->GetFieldEvaluator("overland_coupling")
        ->HasFieldChanged(S_next_.ptr(), "overland_pk");

    Teuchos::RCP <const CompositeVector> source1 =
        S_next_->GetFieldData("overland_source_from_subsurface");

    int c_owned = g->size("cell");
    for (int c=0; c!=c_owned; ++c) {
      (*g)("cell",c) -= (*source1)("cell",c) * (*cell_volume1)("cell",c);
    }
  }
};


} //namespace
} //namespace
