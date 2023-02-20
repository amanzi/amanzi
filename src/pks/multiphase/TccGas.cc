/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Calculates gas component concentration from liquid concentration.
*/

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "TccGas.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
TccGas::TccGas(Teuchos::ParameterList& plist) : MultiphaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  tcc_liquid_key_ = plist_.get<std::string>("tcc liquid key");
  dependencies_.insert(std::make_pair(tcc_liquid_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
TccGas::Clone() const
{
  return Teuchos::rcp(new TccGas(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
TccGas::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& tcc = *S.Get<CompositeVector>(tcc_liquid_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = result_c.MyLength();

  for (int c = 0; c != ncells; ++c) { result_c[0][c] = tcc[0][c] * kH_; }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
TccGas::EvaluatePartialDerivative_(const State& S,
                                   const Key& wrt_key,
                                   const Tag& wrt_tag,
                                   const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = result_c.MyLength();

  if (wrt_key == tcc_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = kH_; }
  }
}

} // namespace Multiphase
} // namespace Amanzi
