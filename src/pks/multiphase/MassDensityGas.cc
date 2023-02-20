/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Calculates mass density of a gas phase.
*/

#include "MassDensityGas.hh"
#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
MassDensityGas::MassDensityGas(Teuchos::ParameterList& plist) : MultiphaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");
  x_gas_key_ = plist_.get<std::string>("mole fraction gas key");
  mol_density_gas_key_ = plist_.get<std::string>("molar density gas key");

  dependencies_.insert(std::make_pair(x_vapor_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(x_gas_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_gas_key_, Tags::DEFAULT));

  mol_mass_H2O_ = plist_.get<double>("molar mass of water");
  mol_mass_ = plist_.get<Teuchos::Array<double>>("molar masses").toVector();
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
MassDensityGas::Clone() const
{
  return Teuchos::rcp(new MassDensityGas(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
MassDensityGas::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  if (S.HasRecord(x_vapor_key_)) {
    const auto& eta_g = *S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");
    auto& result_c = *results[0]->ViewComponent("cell");

    const auto& xv = *S.Get<CompositeVector>(x_vapor_key_).ViewComponent("cell");
    const auto& xg = *S.Get<CompositeVector>(x_gas_key_).ViewComponent("cell");

    int ncells = result_c.MyLength();
    for (int c = 0; c < ncells; ++c) result_c[0][c] = xv[0][c] * mol_mass_H2O_;
  
    for (int i = 0; i < xg.NumVectors(); ++i) {
      for (int c = 0; c < ncells; ++c) {
        result_c[0][c] += xg[i][c] * mol_mass_[i];
      }
    }

    for (int c = 0; c < ncells; ++c) result_c[0][c] *= eta_g[0][c];
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
MassDensityGas::EvaluatePartialDerivative_(const State& S,
                                           const Key& wrt_key,
                                           const Tag& wrt_tag,
                                           const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(false);
}

} // namespace Multiphase
} // namespace Amanzi
