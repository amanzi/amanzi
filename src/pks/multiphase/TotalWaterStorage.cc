/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Multiphase PK

  Field evaluator for total water storage, the conserved quantity:

    TWS = phi * (eta_l * s_l + eta_g * s_g * x_v)
*/

#include "TotalWaterStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalWaterStorage::TotalWaterStorage(Teuchos::ParameterList& plist) : MultiphaseBaseEvaluator(plist)
{
  Init_();
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void
TotalWaterStorage::Init_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  mol_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  mol_density_gas_key_ = plist_.get<std::string>("molar density gas key");
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(saturation_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_gas_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(x_vapor_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalWaterStorage::TotalWaterStorage(const TotalWaterStorage& other)
  : MultiphaseBaseEvaluator(other){};


Teuchos::RCP<Evaluator>
TotalWaterStorage::Clone() const
{
  return Teuchos::rcp(new TotalWaterStorage(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
TotalWaterStorage::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& sl = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& ng = *S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");
  const auto& vg = *S.Get<CompositeVector>(x_vapor_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = phi[0][c] * (sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c] * vg[0][c]);
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
TotalWaterStorage::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& sl = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& ng = *S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");
  const auto& vg = *S.Get<CompositeVector>(x_vapor_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c] * vg[0][c];
    }
  } else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (nl[0][c] - ng[0][c] * vg[0][c]);
    }
  } else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = phi[0][c] * sl[0][c]; }
  } else if (wrt_key == mol_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * vg[0][c]; }
  } else if (wrt_key == x_vapor_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * ng[0][c]; }
  }
}

} // namespace Multiphase
} // namespace Amanzi
