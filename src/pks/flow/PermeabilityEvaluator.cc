/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK
*/

// Amanzi
#include "PDE_HelperDiscretization.hh"

// Amanzi::Flow
#include "FlowDefs.hh"
#include "Permeability_PowerLaw.hh"
#include "Permeability_KozenyCarman.hh"
#include "PermeabilityEvaluator.hh"
#include "PermeabilityModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
PermeabilityEvaluator::PermeabilityEvaluator(Teuchos::ParameterList& plist,
                                             Teuchos::RCP<PermeabilityModelPartition> ppm)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), ppm_(ppm)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(
      std::make_pair(plist_.get<std::string>("permeability porosity factor key"), Tags::DEFAULT));
  }

  // my dependency is pressure.
  porosity_key_ = plist_.get<std::string>("porosity key");
  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
}


PermeabilityEvaluator::PermeabilityEvaluator(const PermeabilityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    ppm_(other.ppm_),
    porosity_key_(other.porosity_key_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
PermeabilityEvaluator::Clone() const
{
  return Teuchos::rcp(new PermeabilityEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PermeabilityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  auto& factor_c = *results[0]->ViewComponent("cell");
  const auto& phi_c = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");

  int ncells = factor_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    factor_c[0][c] = ppm_->second[(*ppm_->first)[c]]->Factor(phi_c[0][c]);
  }

  // optional copy of cell data to boundary
  Operators::CellToBoundaryFaces(*results[0]);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PermeabilityEvaluator::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& results)
{
  auto& factor_c = *results[0]->ViewComponent("cell");
  const auto& phi_c = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");

  int ncells = factor_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    factor_c[0][c] = ppm_->second[(*ppm_->first)[c]]->dFactordPorosity(phi_c[0][c]);
  }

  Operators::CellToBoundaryFaces(*results[0]);
}

} // namespace Flow
} // namespace Amanzi
