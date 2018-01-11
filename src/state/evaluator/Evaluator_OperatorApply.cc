/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! Evaluator_OperatorApply is a secondary evaluator that calculates r = b - Ax

/*!


*/

#include "Evaluator_OperatorApply.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"

namespace Amanzi {

Evaluator_OperatorApply::Evaluator_OperatorApply(Teuchos::ParameterList& plist) :
    EvaluatorSecondary<CompositeVector,CompositeVectorSpace>(plist) {
  // assumes domain-residual_XXX, get the XXX
  Key var_name = Keys::getVarName(my_key_);
  Key operator_name;
  if (var_name.size() >= 10) operator_name = var_name.substr(9,var_name.size());    
  Key domain_name = Keys::getDomain(my_key_);
    
  rhs_key_ = Keys::readKey(plist, domain_name, "rhs", Key("rhs_")+operator_name);
  x_key_ = Keys::readKey(plist, domain_name, "x", Key("x_")+operator_name);

  std::vector<std::string> local_op_keys;
  if (plist.isParameter("local operator keys")) {
    local_op_keys = plist.get<Teuchos::Array<std::string>>("local operator keys").toVector();
  } else {
    auto local_op_suffixes = plist.get<Teuchos::Array<std::string>>("local operator suffixes");
    for (const auto& suffix : local_op_suffixes) {
      local_op_keys.emplace_back(Keys::getKey(domain_name, suffix));
    }
  }

  dependencies_.emplace_back(std::make_pair(rhs_key_, my_tag_));
  dependencies_.emplace_back(std::make_pair(x_key_, my_tag_));
  for (const auto& dep : local_op_keys) {
    dependencies_.emplace_back(std::make_pair(dep, my_tag_));
  }
}

  
// These do the actual work
void
Evaluator_OperatorApply::Evaluate_(const State& S, CompositeVector& result)
{
  const auto& x = S.Get<CompositeVector>(x_key_, my_tag_);
  const auto& b = S.Get<CompositeVector>(rhs_key_, my_tag_);

  // create the global operator
  Operators::Operator_Factory global_op_fac;
  global_op_fac.set_mesh(b.Mesh());
  global_op_fac.set_cvs(x.Map(), b.Map());
  auto global_op = global_op_fac.Create();

  // do the apply
  //  x.Print(std::cout);
  x.ScatterMasterToGhosted();
  result.PutScalarMasterAndGhosted(0.);

  std::cout << "dep lengts = " << dependencies_.size();
  auto op_itr = dependencies_.begin();
  std::cout << "  0th =  " << op_itr->first << std::endl;
  for (std::advance(op_itr, 2); op_itr != dependencies_.end(); ++op_itr) {
    std::cout << "dependency = " << op_itr->first << std::endl;
    S.Get<Operators::Op>(op_itr->first, op_itr->second).ApplyMatrixFreeOp(&*global_op, x, result);
  }

  //  result.Print(std::cout);
  result.Update(1.0, b, -1.0);
}
  
// No derivative
void
Evaluator_OperatorApply::EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Key& wrt_tag, CompositeVector& result)
{
  Errors::Message message;
  message << "Evaluator_OperatorApply for " << my_key_ << " does not implement a derivative.";
  throw(message);
}

} // namespace
