/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "PK_DomainFunctionFactory.hh"

#include "LinearRelaxationEvaluator.hh"

namespace Amanzi {
namespace Evaluators {

LinearRelaxationEvaluator::LinearRelaxationEvaluator(Teuchos::ParameterList& plist,
                                                     Teuchos::RCP<State>& S)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  variable_key_ = plist.get<std::string>("variable key");
  dependencies_.insert(std::make_pair(variable_key_, Tags::DEFAULT));

  auto mesh = S->GetMesh();
  PK_DomainFunctionFactory<PK_DomainFunction> factory(mesh, S);

  auto tmp = plist.sublist("linear relaxation");
  for (auto it = tmp.begin(); it != tmp.end(); ++it) {
    std::string name = it->first;
    if (tmp.isSublist(name)) {
      Teuchos::ParameterList& spec = tmp.sublist(name);
      srcs_.push_back(factory.Create(spec, "source", AmanziMesh::Entity_kind::CELL, Teuchos::null));
    }
  }
}


LinearRelaxationEvaluator::LinearRelaxationEvaluator(const LinearRelaxationEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    variable_key_(other.variable_key_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
LinearRelaxationEvaluator::Clone() const
{
  return Teuchos::rcp(new LinearRelaxationEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
LinearRelaxationEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  double t = S.get_time();
  double dt = S.Get<double>("dt", Tags::DEFAULT);;

  for (int i = 0; i < srcs_.size(); ++i) {
    srcs_[i]->Compute(t, t + dt);
  }

  const auto& var_c = *S.Get<CompositeVector>(variable_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell", false);

  int ncells = result_c.MyLength();
  result_c.PutScalar(0.0);

  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      if (c < ncells) {
        result_c[0][c] = it->second[0] * (var_c[0][c] - it->second[1]);
      }
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
LinearRelaxationEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell", false);

  int ncells = result_c.MyLength();
  result_c.PutScalar(0.0);

  if (wrt_key == variable_key_) {
    double t = S.get_time();
    double dt = S.Get<double>("dt", Tags::DEFAULT);;

    for (int i = 0; i < srcs_.size(); ++i) {
      srcs_[i]->Compute(t, t + dt);
    }

    for (int i = 0; i < srcs_.size(); ++i) {
      for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
        int c = it->first;
        if (c < ncells) result_c[0][c] = it->second[0];
      }
    }
  }
}

}  // namespace Evaluators
}  // namespace Amanzi

