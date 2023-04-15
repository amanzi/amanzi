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

  The porosity model evaluator simply calls the porosity model with
  the correct arguments.
*/

// Amanzi
#include "PDE_HelperDiscretization.hh"

// Amanzi::Flow
#include "FlowDefs.hh"
#include "Porosity_Compressible.hh"
#include "Porosity_Constant.hh"
#include "PorosityEvaluator.hh"
#include "PorosityModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
PorosityEvaluator::PorosityEvaluator(Teuchos::ParameterList& plist,
                                     Teuchos::RCP<PorosityModelPartition> pom)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), pom_(pom)
{
  InitializeFromPlist_();
}


PorosityEvaluator::PorosityEvaluator(const PorosityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pom_(other.pom_),
    pressure_key_(other.pressure_key_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
PorosityEvaluator::Clone() const
{
  return Teuchos::rcp(new PorosityEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void
PorosityEvaluator::InitializeFromPlist_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("porosity key"), Tags::DEFAULT));
  }

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key");
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PorosityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  auto& phi_c = *results[0]->ViewComponent("cell", false);
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->PorosityValue(pres_c[0][c]);
  }

  // optional copy of cell data to boundary
  Operators::CellToBoundaryFaces(*results[0]);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PorosityEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& results)
{
  auto& phi_c = *results[0]->ViewComponent("cell");
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  int ncells = phi_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    phi_c[0][c] = pom_->second[(*pom_->first)[c]]->dPorositydPressure(pres_c[0][c]);
  }

  Operators::CellToBoundaryFaces(*results[0]);
}

} // namespace Flow
} // namespace Amanzi
