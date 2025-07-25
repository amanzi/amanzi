/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Multiphase PK

  Field evaluator for noninear complimentary problem, function G.
*/

#include "NCP_MoleFractions.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
NCP_MoleFractions::NCP_MoleFractions(Teuchos::ParameterList& plist)
  : MultiphaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");
  x_gas_key_ = plist_.get<std::string>("mole fraction gas key");

  dependencies_.insert(std::make_pair(x_vapor_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(x_gas_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_MoleFractions::NCP_MoleFractions(const NCP_MoleFractions& other)
  : MultiphaseEvaluator(other) {};


Teuchos::RCP<Evaluator>
NCP_MoleFractions::Clone() const
{
  return Teuchos::rcp(new NCP_MoleFractions(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
NCP_MoleFractions::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& vg = *S.Get<CompositeVector>(x_vapor_key_).ViewComponent("cell");
  const auto& xg = *S.Get<CompositeVector>(x_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double sum(vg[0][c]);
    for (int i = 0; i < xg.NumVectors() ; ++i) sum += xg[i][c];
    result_c[0][c] = 1.0 - sum;
  }
}


/* ******************************************************************
* Required member: field derivative calculation.
****************************************************************** */
void
NCP_MoleFractions::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == x_vapor_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  } else if (wrt_key == x_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

} // namespace Multiphase
} // namespace Amanzi
