/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for noninear complimentary problem, function G.
*/

#include "NCP_MolarDensities.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
NCP_MolarDensities::NCP_MolarDensities(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");
  molar_density_gas_key_ = plist_.get<std::string>("molar density gas key");
  tcc_gas_key_ = plist_.get<std::string>("tcc gas key");

  dependencies_.insert(std::make_pair(x_vapor_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(molar_density_gas_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(tcc_gas_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_MolarDensities::NCP_MolarDensities(const NCP_MolarDensities& other)
  : MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<Evaluator> NCP_MolarDensities::Clone() const {
  return Teuchos::rcp(new NCP_MolarDensities(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void NCP_MolarDensities::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& xg = *S.Get<CompositeVector>(x_vapor_key_).ViewComponent("cell");
  const auto& ng = *S.Get<CompositeVector>(molar_density_gas_key_).ViewComponent("cell");
  const auto& tcc = *S.Get<CompositeVector>(tcc_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double sum = ng[0][c] * xg[0][c];
    for (int i = 0; i < tcc.NumVectors(); ++i) sum += tcc[i][c];
    result_c[0][c] = ng[0][c] - sum;
  }      
}


/* ******************************************************************
* Required member: field derivative calculation.
****************************************************************** */
void NCP_MolarDensities::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& xg = *S.Get<CompositeVector>(x_vapor_key_).ViewComponent("cell");
  const auto& ng = *S.Get<CompositeVector>(molar_density_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == x_vapor_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -ng[0][c];
  }
  else if (wrt_key == molar_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = 1.0 - xg[0][c];
  }
  else if (wrt_key == tcc_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

