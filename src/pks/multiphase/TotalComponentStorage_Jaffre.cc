/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Multiphase PK

  Field evaluator for a total component storage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (eta_l * s_l + eta_g * s_g)
*/

#include "TotalComponentStorage_Jaffre.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalComponentStorage_Jaffre::TotalComponentStorage_Jaffre(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }

  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  mol_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  mol_density_gas_key_ = plist_.get<std::string>("molar density gas key");

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(saturation_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_gas_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalComponentStorage_Jaffre::TotalComponentStorage_Jaffre(
  const TotalComponentStorage_Jaffre& other)
  : MultiphaseBaseEvaluator(other){};


Teuchos::RCP<Evaluator>
TotalComponentStorage_Jaffre::Clone() const
{
  return Teuchos::rcp(new TotalComponentStorage_Jaffre(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
TotalComponentStorage_Jaffre::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& sl = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& ng = *S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double tmpl = phi[0][c] * nl[0][c] * sl[0][c];
    double tmpg = phi[0][c] * ng[0][c] * (1.0 - sl[0][c]);
    result_c[0][c] = tmpl + tmpg;
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
TotalComponentStorage_Jaffre::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& sl = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& ng = *S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c];
    }
  } else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = phi[0][c] * (nl[0][c] - ng[0][c]); }
  }

  else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = phi[0][c] * sl[0][c]; }
  } else if (wrt_key == mol_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) { result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]); }
  }
}

} // namespace Multiphase
} // namespace Amanzi
