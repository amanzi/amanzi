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
MassDensityGas::MassDensityGas(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  x_gas_key_ = plist_.get<std::string>("mole fraction gas key");
  mol_density_gas_key_ = plist_.get<std::string>("molar density gas key");
  components_ = plist.get<Teuchos::Array<std::string>>("components").toVector();
  mol_mass_ = plist_.get<Teuchos::Array<double>>("molar masses").toVector();

  for (auto& name : components_) {
    Key key = x_gas_key_ + "_" + name;
    dependencies_.insert(std::make_pair(key, Tags::DEFAULT));
  }
  dependencies_.insert(std::make_pair(mol_density_gas_key_, Tags::DEFAULT));

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
  const auto& eta_g = *S.Get<CompositeVector>(mol_density_gas_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  result_c.PutScalar(0.0);

  for (int i = 0; i < components_.size(); ++i) {
    Key key = x_gas_key_ + "_" + components_[i];
    const auto& xg = *S.Get<CompositeVector>(key).ViewComponent("cell");

    for (int c = 0; c < ncells; ++c) {
      result_c[0][c] += xg[0][c] * mol_mass_[i];
    }
  }

  for (int c = 0; c < ncells; ++c) result_c[0][c] *= eta_g[0][c];
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
