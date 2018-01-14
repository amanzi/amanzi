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
#include "Op_Factory.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"

namespace Amanzi {

Evaluator_OperatorApply::Evaluator_OperatorApply(Teuchos::ParameterList& plist) :
    EvaluatorSecondary<CompositeVector,CompositeVectorSpace>(plist) {
  // can this be cleaned up? --etc
  // assumes domain-residual_XXX, get the XXX
  Key var_name = Keys::getVarName(my_key_);
  Key operator_name;
  if (var_name.size() >= 10) operator_name = var_name.substr(9,var_name.size());    
  Key domain_name = Keys::getDomain(my_key_);
    
  rhs_key_ = Keys::readKey(plist, domain_name, "rhs", Key("rhs_")+operator_name);
  x_key_ = Keys::readKey(plist, domain_name, "x", Key("x_")+operator_name);

  if (plist.isParameter("local operator keys")) {
    local_op_keys_ = plist.get<Teuchos::Array<std::string>>("local operator keys").toVector();
  } else {
    auto local_op_suffixes = plist.get<Teuchos::Array<std::string>>("local operator suffixes");
    for (const auto& suffix : local_op_suffixes) {
      local_op_keys_.emplace_back(Keys::getKey(domain_name, suffix));
    }
  }

  // push dependencies
  dependencies_.emplace_back(std::make_pair(rhs_key_, my_tag_));
  dependencies_.emplace_back(std::make_pair(x_key_, my_tag_));
  for (const auto& dep : local_op_keys_) {
    dependencies_.emplace_back(std::make_pair(dep, my_tag_));
  }

  // is domain_map same as range_map?
  is_square_ = plist.get<bool>("square", true);
}

void
Evaluator_OperatorApply::EnsureCompatibility(State& S)
{
  // call base class's to get my requirements set
  EvaluatorSecondary<CompositeVector,CompositeVectorSpace>::EnsureCompatibility(S);

  // set requirements on my dependencies
  //
  // NOTE: this process is a bit wierd.  We know the Mesh of the residual,
  // from the PK, but we don't necessarily know the CVS of the primary
  // variable, since it depends upon the discretization of the Diffusion
  // operator.  This leads us to the following odd setup:
  //
  // 1. PK sets my mesh
  // 2. Push my mesh into the RHS
  // 3. EnsureCompatibility() on the operator -- this sucks in the Mesh from
  //    the RHS, reads the discretization from the PList it gets, and finishes
  //    setting the CVS of the Operator and RHS.
  // 4. Now we know the schema of the RHS, use it to set the RHS of the
  //    primary variable and the residual.
  auto& dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_);
  if (dep_fac.Mesh() != Teuchos::null) {
    // -- set the rhs's mesh
    {
      auto& rhs_fac = S.Require<CompositeVector,CompositeVectorSpace>(rhs_key_, my_tag_);
      rhs_fac.Update(dep_fac);
    }
    
    // -- call ensure compatibility of all dependencies, propagating info, and in
    //    the case of local ops and rhs, setting their schema
    for (auto& dep : dependencies_) {
      S.RequireEvaluator(dep.first, dep.second)->EnsureCompatibility(S);
    }

    // -- re-get the rhs face, which now includes schema info
    // copy intentional
    auto rhs_fac = S.Require<CompositeVector,CompositeVectorSpace>(rhs_key_, my_tag_);
    rhs_fac.SetOwned(false);
    
    if (is_square_) {
      // -- x schema is the same map as mine
      auto& x_fac = S.Require<CompositeVector,CompositeVectorSpace>(x_key_, my_tag_);
      x_fac.Update(rhs_fac);
    }

    // update my schema
    S.Require<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_, my_key_)
        .Update(rhs_fac);
    
    // -- Now we can check the domain and range
    // -- range space of operator is a subset of mine
    // THIS SHOULD BE UNNECESSARY, as we got the schema from diffusion.  Idiot check!
    for (const auto& op : local_op_keys_) {
      auto& op_fac = S.Require<Operators::Op, Operators::Op_Factory>(op, my_tag_);
      if (!Operators::cvsFromSchema(op_fac.schema_col(), op_fac.mesh()).SubsetOf(rhs_fac)) {
        Errors::Message message;
        message << "Evaluator_OperatorApply for " << my_key_ << ": Op " << op << " column schema is not a subset of my schema.";
        throw(message);
      }

      if (is_square_) {
        if (!Operators::cvsFromSchema(op_fac.schema_row(), op_fac.mesh()).SubsetOf(rhs_fac)) {
          Errors::Message message;
          message << "Evaluator_OperatorApply for " << my_key_ << ": Op " << op << " row schema is not a subset of my schema.";
          throw(message);
        }
      }
    }
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
  x.ScatterMasterToGhosted();
  result.PutScalarMasterAndGhosted(0.);

  auto op_itr = dependencies_.begin();
  for (std::advance(op_itr, 2); op_itr != dependencies_.end(); ++op_itr) {
    S.Get<Operators::Op>(op_itr->first, op_itr->second).ApplyMatrixFreeOp(&*global_op, x, result);
  }

  // b - Ax
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
