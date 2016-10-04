/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "snow_distribution.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void SnowDistribution::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(S_next_.ptr());

  // update the stiffness matrix
  matrix_->Init();
  Teuchos::RCP<const CompositeVector> cond =
<<<<<<< HEAD
    S_next_->GetFieldData(getKey(domain_,"upwind_snow_conductivity"), name_);
=======
    S_next_->GetFieldData("upwind_snow_conductivity", name_);
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  matrix_diff_->ApplyBCs(true, true);
  
  // update the potential
<<<<<<< HEAD
  S->GetFieldEvaluator(getKey(domain_,"snow_skin_potential"))->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> potential = S->GetFieldData(getKey(domain_,"snow_skin_potential"));
=======
  S->GetFieldEvaluator("snow_skin_potential")->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> potential = S->GetFieldData("snow_skin_potential");
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  // calculate the residual
  matrix_->ComputeNegativeResidual(*potential, *g);
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void SnowDistribution::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  //  double dt = S_next_->time() - S_inter_->time();

  // get these fields
  Teuchos::RCP<const CompositeVector> h1 = S_next_->GetFieldData(key_);
  Teuchos::RCP<const CompositeVector> cv1 =
<<<<<<< HEAD
    S_next_->GetFieldData("surface-cell_volume");
=======
    S_next_->GetFieldData("surface_cell_volume");
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  std::vector<double> time(1,S_next_->time());
  double precip = (*precip_func_)(time);
  
  g->ViewComponent("cell",false)->Multiply(10.,
          *cv1->ViewComponent("cell",false), *h1->ViewComponent("cell",false), 1.);
  g->ViewComponent("cell",false)->Update(-precip*10., *cv1->ViewComponent("cell",false), 1.);

};


} //namespace
} //namespace
