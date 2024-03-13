/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Evaluator for volumetric flow rate.
*/

#include "CommonDefs.hh"
#include "Key.hh"
#include "PDE_HelperDiscretization.hh"
#include "Upwind.hh"

#include "VolumetricFlowRateEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor
****************************************************************** */
VolumetricFlowRateEvaluator::VolumetricFlowRateEvaluator(Teuchos::ParameterList& plist,
                                                         double molar_rho)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), molar_rho_(molar_rho)
{
  AMANZI_ASSERT(my_keys_.size() > 0);

  domain_ = Keys::getDomain(my_keys_[0].first);
  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");
  dependencies_.insert(std::make_pair(mol_flowrate_key_, Tags::DEFAULT));

  if (molar_rho == 0.0) {
    mol_density_key_ = Keys::getKey(domain_, "molar_density_liquid");
    dependencies_.insert(std::make_pair(mol_density_key_, Tags::DEFAULT));
  }
}


/* ******************************************************************
* A copy constructor.
****************************************************************** */
VolumetricFlowRateEvaluator::VolumetricFlowRateEvaluator(const VolumetricFlowRateEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    mol_flowrate_key_(other.mol_flowrate_key_),
    upwind_(other.upwind_),
    bc_(other.bc_)
{}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<Evaluator>
VolumetricFlowRateEvaluator::Clone() const
{
  return Teuchos::rcp(new VolumetricFlowRateEvaluator(*this));
}


/* ******************************************************************
* Evaluator
****************************************************************** */
void
VolumetricFlowRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // BCs may not be available during initialization time
  if (bc_.get()) {
    std::vector<int>& bc_model = bc_->bc_model();
    const auto& mol_flowrate = S.Get<CompositeVector>(mol_flowrate_key_);
    const auto& mol_flowrate_f = *mol_flowrate.ViewComponent("face");

    auto& result_f = *results[0]->ViewComponent("face");
    int nfaces = result_f.MyLength();

    if (upwind_.get()) {
      const auto& mol_density = S.Get<CompositeVector>(mol_density_key_);

      // populate auxiliary field for upwind purpose
      CompositeVectorSpace cvs = mol_flowrate.Map();
      cvs.SetOwned(false);
      cvs.AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

      CompositeVector beta(cvs);
      *beta.ViewComponent("cell") = *mol_density.ViewComponent("cell");
      if (mol_density.HasComponent("boundary_face")) {
        Operators::BoundaryFacesToFaces(bc_model, mol_density, beta);
      } else {
        Operators::CellToBoundaryFaces(bc_model, beta);
      }

      upwind_->Compute(mol_flowrate, bc_model, beta);

      auto& beta_f = *beta.ViewComponent("face");
      for (int f = 0; f < nfaces; ++f) result_f[0][f] = mol_flowrate_f[0][f] / beta_f[0][f];
    } else {
      for (int f = 0; f < nfaces; ++f) result_f[0][f] = mol_flowrate_f[0][f] / molar_rho_;
    }
  } else {
    results[0]->PutScalar(0.0);
  }
}

} // namespace Flow
} // namespace Amanzi
