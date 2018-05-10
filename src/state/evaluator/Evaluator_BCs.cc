/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! Evaluates a function to provide Dirichlet data on faces.

#include "Evaluator_BCs.hh"

namespace Amanzi {

Evaluator_BCs::Evaluator_BCs(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<Operators::BCs, Operators::BCs_Factory>(plist) {
  for (auto sublist : plist.sublist("boundary functions")) {
    std::string sublist_name = Keys::cleanPListName(sublist.first);
    dependencies_.push_back(std::make_pair(sublist_name, my_tag_));
    std::string bc_type = plist.sublist("boundary functions")
                          .sublist(sublist_name).get<std::string>("boundary condition type");
    if (bc_type == "Dirichlet")
      bc_types_.push_back(Operators::OPERATOR_BC_DIRICHLET);
    else if (bc_type == "Neumann")
      bc_types_.push_back(Operators::OPERATOR_BC_NEUMANN);
    else if (bc_type == "mixed")
      bc_types_.push_back(Operators::OPERATOR_BC_MIXED);
    else {
      Errors::Message msg;
      msg << "BC for " << my_key_ << " has unknown type \"" << bc_type << "\"";
      throw(msg);
    }
  }
}

bool Evaluator_BCs::IsDifferentiableWRT(const State& S, const Key& wrt_key,
                                        const Key& wrt_tag) const {
  return false;
}

// void Evaluator_BCs::EnsureCompatibleDerivative(State &S,
//         const Key &wrt_key, const Key &wrt_tag) {
//   Errors::Message msg("BCs are not differentiable");
//   throw(msg);
// }

void Evaluator_BCs::EnsureCompatibility(State& S) {
  auto& my_fac = S.Require<Operators::BCs, Operators::BCs_Factory>(my_key_, my_tag_, my_key_);
  if (my_fac.mesh().get()) {
    for (const auto& dep : dependencies_) {
      auto& eval = S.RequireEvaluator(dep.first, dep.second);
      auto& fac = S.Require<Functions::BoundaryFunction, Functions::BoundaryFunctionFactory>(
          dep.first, dep.second);
      fac.set_mesh(my_fac.mesh());
      auto& bf_list = plist_.sublist("boundary functions").sublist(dep.first);
      fac.set_parameterlist(bf_list);
      eval.EnsureCompatibility(S);      
    }
  }
}

void Evaluator_BCs::Evaluate_(const State &S, Operators::BCs&result) {
  auto& model = result.bc_model();
  auto& value = result.bc_value();

  int i = 0;
  for (const auto& dep : dependencies_) {
    const auto& bc_value = S.Get<Functions::BoundaryFunction>(dep.first, dep.second);
    for (const auto& fv : bc_value) {
      model[fv.first] = bc_types_[i];
      value[fv.first] = fv.second;
    }
  }
}

}  // namespace Amanzi
