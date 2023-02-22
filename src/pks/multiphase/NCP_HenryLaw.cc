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

#include "NCP_HenryLaw.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
NCP_HenryLaw::NCP_HenryLaw(Teuchos::ParameterList& plist) : MultiphaseEvaluator(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  pressure_gas_key_ = plist_.get<std::string>("pressure gas key");
  mol_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");

  dependencies_.insert(std::make_pair(pressure_gas_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_HenryLaw::NCP_HenryLaw(const NCP_HenryLaw& other) : MultiphaseEvaluator(other){};


Teuchos::RCP<Evaluator>
NCP_HenryLaw::Clone() const
{
  return Teuchos::rcp(new NCP_HenryLaw(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
NCP_HenryLaw::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& pg = *S.Get<CompositeVector>(pressure_gas_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  for (int c = 0; c != ncells; ++c) { result_c[0][c] = pg[0][c] * kH_ - nl[0][c]; }
}


/* ******************************************************************
* Required member: field derivative calculation.
****************************************************************** */
void
NCP_HenryLaw::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  if (wrt_key == pressure_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = kH_;
  } else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

} // namespace Multiphase
} // namespace Amanzi
