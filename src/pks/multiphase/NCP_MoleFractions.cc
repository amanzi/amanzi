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
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }

  x_key_ = plist_.get<std::string>("mole fraction key");
  components_ = plist.get<Teuchos::Array<std::string>>("components").toVector();

  for (auto& name : components_) {
    Key key = x_key_ + "_" + name;
    dependencies_.insert(std::make_pair(key, Tags::DEFAULT));
  }
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_MoleFractions::NCP_MoleFractions(const NCP_MoleFractions& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other) {};


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
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell", false);

  result_c.PutScalar(1.0);
  for (auto& name : components_) {
    Key key = x_key_ + "_" + name;
    const auto& xg = *S.Get<CompositeVector>(key).ViewComponent("cell");

    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] -= xg[0][c];
    }
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

  for (auto& name : components_) {
    Key key = x_key_ + "_" + name;
    if (wrt_key == key) {
      for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
    }
  }
}

} // namespace Multiphase
} // namespace Amanzi
