/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Field evaluator for conserved quantity: the total component storage

    TCS = phi * (s_l * C_l + s_g * C_g)

  where C_p is the concentration a component in phase p.
*/

#include "TotalComponentStorage_Tcc.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalComponentStorage_Tcc::TotalComponentStorage_Tcc(Teuchos::ParameterList& plist) :
    MultiphaseBaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }

  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  tcc_liquid_key_ = plist_.get<std::string>("total component concentration liquid key");
  tcc_gas_key_ = plist_.get<std::string>("total component concentration gas key");

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(saturation_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(tcc_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(tcc_gas_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalComponentStorage_Tcc::TotalComponentStorage_Tcc(const TotalComponentStorage_Tcc& other) :
    MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<Evaluator> TotalComponentStorage_Tcc::Clone() const {
  return Teuchos::rcp(new TotalComponentStorage_Tcc(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorage_Tcc::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& sl = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& tccl = *S.Get<CompositeVector>(tcc_liquid_key_).ViewComponent("cell");
  const auto& tccg = *S.Get<CompositeVector>(tcc_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = phi[0][c] * (tccl[0][c] * sl[0][c] + tccg[0][c] * (1.0 - sl[0][c]));
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorage_Tcc::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& sl = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& tccl = *S.Get<CompositeVector>(tcc_liquid_key_).ViewComponent("cell");
  const auto& tccg = *S.Get<CompositeVector>(tcc_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = tccl[0][c] * sl[0][c] + tccg[0][c] * (1.0 - sl[0][c]);
    }
  }
  else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (tccl[0][c] - tccg[0][c]);
    }
  }

  else if (wrt_key == tcc_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c];
    }

  } else if (wrt_key == tcc_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]);
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
