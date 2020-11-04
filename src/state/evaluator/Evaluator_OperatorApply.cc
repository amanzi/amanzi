/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Evaluator_OperatorApply is a secondary evaluator that calculates r = b - Ax

/*!

This is expected to be the residual equation for a Leaf PK, therefore it lives
on a domain.

It is expected to be of the form:

r = rhs_1 + ... rhs_k + rhs_A00 + rhs_A01 + ... rhs_Anm - (A00 x0 + A01 x0 + ...
A0n x0 + A11 x1 + ... Anm xn)

Where:

  x0                   | is the primary variable associated with this residual
  x1 ... xn            | are off-diagonal primary variables.
  A0i, rhs_A0i         | are diagonal, local operators and their RHS boundary
                       | conditions/sources.
  Aji, rhs_Aji         | are not-necessarily diagonal, local operators and their
                       | RHSs applied to the xj.
  rhs_1 ... rhs_k      | are arbitrary source terms which CANNOT NOT BE AFFECTED
BY | boundary conditions (i.e, for conservation equations | discretized using
control volume methods, BCs affect only | faces while sources are on cells.

Note that we can infer some constraints here:

- The domain and range of the A0i must be subsets of the r space.
  Realistically, r's space is the union of the domain and range spaces of the
  A0i and the space of the rhs_k.

-


*/

#include "Evaluator_OperatorApply.hh"
#include "Op_Factory.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"

namespace Amanzi {

Evaluator_OperatorApply::Evaluator_OperatorApply(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  Key domain = Keys::getDomain(my_keys_[0].first);

  // x, x0
  x0_key_ = Keys::readKey(plist_, domain, "diagonal primary x");
  Teuchos::Array<Key> empty;
  x_keys_ = Keys::readKeys(plist_, domain, "off-diagonal primary x", &empty);

  // op0, op0 rhs's
  op0_keys_ =
    Keys::readKeys(plist_, domain, "diagonal local operators", &empty);
  {
    Teuchos::Array<Key> defaults(op0_keys_.size());
    int i = 0;
    for (const auto& op_key : op0_keys_) {
      defaults[i] = op_key + "_rhs";
      ++i;
    }
    op0_rhs_keys_ =
      Keys::readKeys(plist_, domain, "diagonal local operator rhss", &defaults);
  }

  // opi, opi rhs's
  int j = 0;
  for (const auto& x_key : x_keys_) {
    // strip the domain from x_key
    Key x_key_suffix = Keys::getVarName(x_key);
    op_keys_.emplace_back(Keys::readKeys(
      plist_, domain, x_key_suffix + " local operators"));

    Teuchos::Array<Key> defaults(op_keys_[j].size());
    int i = 0;
    for (const auto& op_key : op_keys_[j]) {
      defaults[i] = op_key + "_rhs";
      ++i;
    }
    op_rhs_keys_.emplace_back(Keys::readKeys(
      plist_, domain, x_key_suffix + " local operator rhss", &defaults));
  }

  // other rhss
  rhs_keys_ = Keys::readKeys(plist_, domain, "additional rhss", &empty);
  if (rhs_keys_.size() > 0) {
    rhs_scalars_ = plist_.get<Teuchos::Array<double>>("rhs coefficients");
    if (rhs_keys_.size() != rhs_scalars_.size()) {
      Errors::Message message;
      message << "Evaluator_OperatorApply for " << my_keys_[0].first
              << ": RHS terms names (" << (int)rhs_keys_.size()
              << ") and scalings (" << (int)rhs_scalars_.size()
              << ") are not the same length";
      throw(message);
    }
  }

  // get a key to store the global operator
  Key x0_key_suffix = Keys::getVarName(x0_key_);
  global_operator_key_ = Keys::readKey(plist_, domain,
          x0_key_suffix + " global operator",
          x0_key_suffix + "_global_operator");
  
  // push dependencies
  dependencies_.emplace_back(std::make_pair(x0_key_, my_keys_[0].second));
  for (const auto& x_key : x_keys_)
    dependencies_.emplace_back(std::make_pair(x_key, my_keys_[0].second));
  for (const auto& rhs : rhs_keys_)
    dependencies_.emplace_back(std::make_pair(rhs, my_keys_[0].second));
  for (const auto& op_key : op0_keys_)
    dependencies_.emplace_back(std::make_pair(op_key, my_keys_[0].second));
  for (const auto& op_rhs_key : op0_rhs_keys_)
    dependencies_.emplace_back(std::make_pair(op_rhs_key, my_keys_[0].second));
  for (const auto& op_list : op_keys_)
    for (const auto& op_key : op_list)
      dependencies_.emplace_back(std::make_pair(op_key, my_keys_[0].second));
  for (const auto& rhs_list : op_rhs_keys_)
    for (const auto& rhs_key : rhs_list)
      dependencies_.emplace_back(std::make_pair(rhs_key, my_keys_[0].second));

  // primary entity is the entity on which the operator primarily exists.
  // This is CELL for control volume methods and NODE for Lagrangian FE
  // methods.
  //
  // All of x, residual, and rhs must have this entity.
  primary_entity_ = plist_.get<std::string>("primary entity", "cell");
  primary_entity_kind_ = AmanziMesh::entity_kind(primary_entity_);
}

void
Evaluator_OperatorApply::EnsureCompatibility(State& S)
{
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_keys_.size() == 1);
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(
      my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);
  for (const auto& dep : dependencies_)
    S.RequireEvaluator(dep.first, dep.second);

  // check plist for vis or checkpointing control
  EnsureCompatibility_Flags_(S);

  // set requirements on my dependencies
  //
  // NOTE: this process is a bit wierd.  We know the Mesh of the residual,
  // from the PK, but we don't necessarily know the CVS of the primary
  // variable, since it depends upon the discretization of the operator.  This
  // leads us to the following odd setup:
  //
  // 1. PK sets my mesh
  // 2. Push my mesh into all local operators' RHSs
  // 3. EnsureCompatibility() on the operator/rhs -- this sucks in the Mesh from
  //    the RHS, reads the discretization from the PList it gets, and finishes
  //    setting the CVS of the Operator and RHS.
  // 4. Now we know the space of the RHSs, use it to set the space of the
  //    primary variable and the residual.
  //
  // Futhermore, derivatives must be considered.
  //
  bool has_derivs = S.HasDerivativeSet(my_keys_[0].first, my_keys_[0].second);
  if (my_fac.Mesh() != Teuchos::null && my_fac.size() == 0) {
    db_ = Teuchos::rcp(new Debugger(my_fac.Mesh(), std::string("OperatorApply: ")+my_keys_[0].first, plist_));

    // add the primary
    my_fac.AddComponent(primary_entity_, primary_entity_kind_, 1);
    CompositeVectorSpace my_fac_mesh_only(my_fac);

    // -- set the rhss' mesh, need for primary entity
    for (const auto& rhs_key : op0_rhs_keys_) {
      auto& rhs_fac = S.Require<CompositeVector, CompositeVectorSpace>(
        rhs_key, my_keys_[0].second);
      rhs_fac.Update(my_fac_mesh_only);
    }

    for (const auto& rhs_list : op_rhs_keys_) {
      for (const auto& rhs_key : rhs_list) {
        auto& rhs_fac = S.Require<CompositeVector, CompositeVectorSpace>(
          rhs_key, my_keys_[0].second);
        rhs_fac.Update(my_fac_mesh_only);
      }
    }

    // -- call ensure compatibility of all op dependencies, setting their mesh
    //    and getting their schema
    for (const auto& rhs_key : op0_rhs_keys_)
      S.RequireEvaluator(rhs_key, my_keys_[0].second).EnsureCompatibility(S);
    for (const auto& rhs_list : op_rhs_keys_)
      for (const auto& rhs_key : rhs_list)
        S.RequireEvaluator(rhs_key, my_keys_[0].second).EnsureCompatibility(S);

    // -- determine the space of the residual as the union of spaces of the
    //    op's rhs
    for (const auto& rhs_key : op0_rhs_keys_)
      my_fac.Update(S.Require<CompositeVector, CompositeVectorSpace>(
        rhs_key, my_keys_[0].second));
    for (const auto& rhs_list : op_rhs_keys_)
      for (const auto& rhs_key : rhs_list)
        my_fac.Update(S.Require<CompositeVector, CompositeVectorSpace>(
          rhs_key, my_keys_[0].second));

    // -- x must be the same as r
    S.Require<CompositeVector, CompositeVectorSpace>(x0_key_,
                                                     my_keys_[0].second)
      .Update(my_fac);
    S.RequireEvaluator(x0_key_, my_keys_[0].second).EnsureCompatibility(S);

    // -- require the global op for forward application/Update_ calls
    auto& global_op_fac = S.Require<Operators::Operator, Operators::Operator_Factory>(
        global_operator_key_, my_keys_[0].second, global_operator_key_);
    global_op_fac.set_mesh(my_fac.Mesh());
    global_op_fac.set_cvs(my_fac.CreateSpace());
    global_op_fac.set_plist(Teuchos::rcp(new Teuchos::ParameterList(plist_)));

    // NOTE: should we ensure that other x's are a superset of the
    // domain-space of the other ops? Otherwise there is nothing to say about
    // the other x's.  --etc

    // -- other RHSs must be a subset of r's space
    for (const auto& rhs_key : rhs_keys_) {
      S.Require<CompositeVector, CompositeVectorSpace>(rhs_key,
                                                       my_keys_[0].second)
        .Update(my_fac_mesh_only);
      S.RequireEvaluator(rhs_key, my_keys_[0].second).EnsureCompatibility(S);
      if (!S.Require<CompositeVector, CompositeVectorSpace>(rhs_key,
                                                            my_keys_[0].second)
             .SubsetOf(my_fac)) {
        Errors::Message message;
        message << "Evaluator_OperatorApply for " << my_keys_[0].first
                << ": RHS " << rhs_key << " space is not a subset of my space.";
        throw(message);
      }
    }

    if (has_derivs) {
      for (const auto& deriv :
           S.GetDerivativeSet(my_keys_[0].first, my_keys_[0].second)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        AMANZI_ASSERT(wrt.second == my_keys_[0].second);
        AMANZI_ASSERT(
          wrt.first ==
          x0_key_); // NEED TO IMPLEMENT OFF-DIAGONALS EVENTUALLY --etc
        auto& deriv_fac =
          S.RequireDerivative<Operators::Operator, Operators::Operator_Factory>(
            my_keys_[0].first,
            my_keys_[0].second,
            wrt.first,
            wrt.second,
            my_keys_[0].first);
        deriv_fac.set_mesh(my_fac.Mesh());

        // FIX ME TO BE NON_SQUARE FOR OFF DIAGONALS? --etc
        deriv_fac.set_cvs(my_fac.CreateSpace());
        deriv_fac.set_plist(Teuchos::rcp(new Teuchos::ParameterList(plist_)));
      }
    }
  }

  // require data for derivatives
  if (has_derivs) {
    for (const auto& deriv :
         S.GetDerivativeSet(my_keys_[0].first, my_keys_[0].second)) {
      auto wrt = Keys::splitKeyTag(deriv.first);
      AMANZI_ASSERT(wrt.second == my_keys_[0].second);
      AMANZI_ASSERT(wrt.first == x0_key_); // NEED TO IMPLEMENT OFF-DIAGONALS EVENTUALLY --etc

      // quasi-linear operators covered by Update() call, no need to require anything

      // jacobian terms
      for (const auto& op_key : op0_keys_) {
        if (plist_.get<bool>("include Jacobian term d" + op_key + "_d" + wrt.first, false)) {
          // there is a Jacobian term available
          S.RequireDerivative<Operators::Op, Operators::Op_Factory>(
              op_key, my_keys_[0].second, wrt.first, wrt.second);
          S.RequireEvaluator(op_key, my_keys_[0].second)
              .EnsureCompatibility(S);
        }
      }

      // rhs terms
      for (const auto& rhs_key : rhs_keys_) {
        if (S.GetEvaluator(rhs_key, my_keys_[0].second)
              .IsDifferentiableWRT(S, wrt.first, wrt.second)) {
          // require the derivative, and update its factory to match its value
          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
             rhs_key, my_keys_[0].second, wrt.first, wrt.second)
            .Update(S.Require<CompositeVector, CompositeVectorSpace>(
              rhs_key, my_keys_[0].second));
          // re-call EnsureCompatibility() to make sure these derivative
          // requirements get processed
          S.RequireEvaluator(rhs_key, my_keys_[0].second)
            .EnsureCompatibility(S);
        }
      }
    }
  }
}



// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
bool
Evaluator_OperatorApply::UpdateDerivative(State& S, const Key& requestor,
        const Key& wrt_key, const Key& wrt_tag)
{
  if (!IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
    Errors::Message msg;
    msg << "EvaluatorSecondary (" << my_keys_[0].first << "," << my_keys_[0].second << ") is not differentiable with respect to (" << wrt_key << "," << wrt_tag << ").";
    throw(msg);
  }

  Key wrt = Keys::getKeyTag(wrt_key, wrt_tag);
  auto& deriv_request_set = deriv_requests_[wrt];
  
  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Evaluator_OperatorApply Derivative d" << my_keys_[0].first << "_d" << wrt_key
              << " requested by " << requestor << "..." << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = deriv_request_set.empty(); // not done once

  // -- must update if our our dependencies have changed, as these affect the
  // partial derivatives
  Key my_request = Key{ "d" } +
                   Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second) +
                   "_d" + Keys::getKeyTag(wrt_key, wrt_tag);

  // I don't believe this should be required.  It may be a lot of extra work to
  //  compute ourselves (e.g. what if this is assembling local matrices).
  //  update |= Update(S, my_request);

  // -- must update if any of our dependencies or our dependencies' derivatives have changed
  // check x0 key
  Key my_tag = my_keys_[0].second;
  update |= S.GetEvaluator(x0_key_, my_tag).Update(S, my_request);

  // check rhs keys
  for (const auto& dep : rhs_keys_) {
    // but we must update our dependencies as we will use them to compute our derivative
    update |= S.GetEvaluator(dep, my_tag).Update(S, my_request);

    // and we must update our dependencies derivatives to apply the chain rule
    if (S.GetEvaluator(dep, my_tag)
          .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
      update |= S.GetEvaluator(dep, my_tag)
                  .UpdateDerivative(S, my_request, wrt_key, wrt_tag);
    }
  }

  // check op keys
  for (const auto& dep : op0_keys_) {
    update |= S.GetEvaluator(dep, my_tag).Update(S, my_request); // updates the primary term
    if (plist_.get<bool>("include Jacobian term d" + dep + "_d" + wrt_key, false)) {
      // check to see if we include the Jacobian term
      if (S.GetEvaluator(dep, my_tag)
          .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
        update |= S.GetEvaluator(dep, my_tag)
                  .UpdateDerivative(S, my_request, wrt_key, wrt_tag);
      }
    }
  }

  // placeholder for op -- off-diagonal keys
  // check offdiagonal x keys
  for (const auto& dep : x_keys_) {
    // but we must update our dependencies as we will use them to compute our derivative
    update |= S.GetEvaluator(dep, my_tag).Update(S, my_request);

    // and we must update our dependencies derivatives to apply the chain rule
    if (S.GetEvaluator(dep, my_tag)
          .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
      update |= S.GetEvaluator(dep, my_tag)
                  .UpdateDerivative(S, my_request, wrt_key, wrt_tag);
    }
  }
  // check offdiagonal op keys
  for (const auto& op_list : op_keys_) {
    for (const auto& dep : op_list) {
      update |= S.GetEvaluator(dep, my_tag).Update(S, my_request); // updates the primary term
      if (plist_.get<bool>("include Jacobian term d" + dep + "_d" + wrt_key, false)) {
        // check to see if we include the Jacobian term
        if (S.GetEvaluator(dep, my_tag)
            .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
          update |= S.GetEvaluator(dep, my_tag)
                    .UpdateDerivative(S, my_request, wrt_key, wrt_tag);
        }
      }
    }
  }

  // Do the update
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "... updating." << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateDerivative_(S, wrt_key, wrt_tag);
    deriv_request_set.clear();
    deriv_request_set.insert(requestor);
    return true;
  } else {
    // Otherwise, simply service the request
    if (deriv_request_set.find(requestor) == deriv_request_set.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "... not updating but new to this request." << std::endl;
      }
      deriv_request_set.insert(requestor);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "... has not changed." << std::endl;
      }
      return false;
    }
  }
}



// These do the actual work
void
Evaluator_OperatorApply::Update_(State& S)
{
  Key my_tag = my_keys_[0].second;
  auto& result = S.GetW<CompositeVector>(
    my_keys_[0].first, my_tag, my_keys_[0].first);

  const auto& x = S.Get<CompositeVector>(x0_key_, my_tag);
  x.ScatterMasterToGhosted();
  db_->WriteVector("x", x);
  result.putScalarMasterAndGhosted(0.);

  if (S.Get<Operators::Operator>(global_operator_key_, my_tag).size() == 0) {
    // first pass through -- set up the global operator by giving it all the
    // local ops
    auto& global_operator = S.GetW<Operators::Operator>(global_operator_key_,
            my_tag, global_operator_key_);
    for (const auto& op_key : op0_keys_) {
      // fix const correctness sometime?  Likely this would require a redesign...
      global_operator.OpPushBack(S.GetPtrW<Operators::Op>(op_key, my_tag, op_key));
    }
  }    

  // This is still a bit ugly -- each op has its own RHS vector, because BCs
  // were not applied globally.  So we have to call the forward apply, then
  // deal with RHSs independently
  //
  // Fixing this requires a rewrite of Operator to make it more functional and
  // split out Op and Operator. --etc
  const auto& global_operator = S.Get<Operators::Operator>(global_operator_key_, my_tag);
  AMANZI_ASSERT(global_operator.size() == op0_keys_.size());
  global_operator.apply(x, result);
  db_->WriteVector(std::string("A * ")+x0_key_, result);
  
  for (const auto& op_rhs_key : op0_rhs_keys_) {
    // Ax - b
    result.update(-1.0, S.Get<CompositeVector>(op_rhs_key, my_tag), 1.0);
  }
  db_->WriteVector(std::string("A * ")+x0_key_+" -b", result);

  int j = 0;
  for (const auto& op_list : op_keys_) {
    const auto& xj = S.Get<CompositeVector>(x_keys_[j], my_tag);
    xj.ScatterMasterToGhosted();

    int k = 0;
    for (const auto& op_key : op_list) {
      // FIX ME BEFORE USING -- store a Global operator for each in the same
      // way we did for x0 --etc
      
      // create the global operator
      Operators::Operator_Factory global_op_fac;
      global_op_fac.set_mesh(xj.Mesh());
      global_op_fac.set_cvs(xj.getMap(), result.getMap());
      auto global_op = global_op_fac.Create();

      // do the apply
      S.Get<Operators::Op>(op_key, my_tag)
        .ApplyMatrixFreeOp(&*global_op, xj, result);
      db_->WriteVector(op_key+" Ax", result);

      // Ax - b
      result.update(
        -1.0,
        S.Get<CompositeVector>(op_rhs_keys_[j][k], my_tag),
        1.0);
      db_->WriteVector(op_key+" Ax-b", result);
      ++k;
    }
    ++j;
  }

  // add all the additional RHSs
  j = 0;
  for (const auto& rhs_key : rhs_keys_) {
    result.update(rhs_scalars_[j],
                  S.Get<CompositeVector>(rhs_key, my_tag),
                  1.0);
    db_->WriteVector(rhs_key, result);
    ++j;
  }
}


// creates the operator for applying inverses
//
// Things to note:
//  - Update_() has been called, so we have the local operators.
//  - UpdateDerivative(dep, wrt) has been called on all dependencies that depend
//  upon wrt.
//
void
Evaluator_OperatorApply::UpdateDerivative_(State& S, const Key& wrt_key,
                                           const Key& wrt_tag)
{
  AMANZI_ASSERT(wrt_tag == my_keys_[0].second);
  auto global_op = S.GetDerivativePtrW<Operators::Operator>(
    my_keys_[0].first, my_keys_[0].second, wrt_key, wrt_tag, my_keys_[0].first);

  if (global_op->size() == 0) {
    // push in local ops the first time
    if (wrt_key == x0_key_) {
      // diagonal entry
      // collect all operators and jacobian info
      for (const auto& op_key : op0_keys_) {
        // FIX ME AND MAKE THIS CONST CORRECT --etc
        //std::cout << "Adding diffusion op to operator" << std::endl;
        global_op->OpPushBack(
          S.GetPtrW<Operators::Op>(op_key, my_keys_[0].second, op_key));
        if (S.GetEvaluator(op_key, my_keys_[0].second).IsDifferentiableWRT(S, wrt_key, wrt_tag) &&
            plist_.get<bool>("include Jacobian term d" + op_key + "_d" + wrt_key, false)) {
          global_op->OpPushBack(S.GetDerivativePtrW<Operators::Op>(
            op_key, my_keys_[0].second, wrt_key, wrt_tag, op_key));
        }
      }

      int j = 0;
      for (const auto& rhs_key : rhs_keys_) {
        if (S.GetEvaluator(rhs_key, my_keys_[0].second)
              .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
          // FIX ME AND MAKE THIS LESS HACKY AND CONST CORRECT --etc
          CompositeVector drhs = S.GetDerivativeW<CompositeVector>(
            rhs_key, my_keys_[0].second, wrt_key, wrt_tag, rhs_key);
          for (const auto& comp : drhs) {
            if (AmanziMesh::entity_kind(comp) == AmanziMesh::CELL) {
              auto op_cell =
                Teuchos::rcp(new Operators::Op_Cell_Cell(rhs_key, drhs.Mesh()));
              global_op->OpPushBack(op_cell);
            } else {
              AMANZI_ASSERT(0);
            }
          }
        }
        ++j;
      }
    } else {
      // off diagonal
      throw("OperatorApply not implemented yet for offdiagonals!");
    }

    // symbolic assemble the first time
    //    global_op->SymbolicAssembleMatrix();
    auto plist = Teuchos::parameterList(plist_.sublist("preconditioner"));
    global_op->set_inverse_parameters(*plist);
  }

  // push values
  if (wrt_key == x0_key_) {
    // count following the above OpPushBack calls to find my Op
    int i = 0;
    for (const auto& op_key : op0_keys_) {
      i++;
      if (S.GetEvaluator(op_key, my_keys_[0].second).IsDifferentiableWRT(S, wrt_key, wrt_tag) &&
          plist_.get<bool>("include Jacobian term d" + op_key + "_d" + wrt_key, false)) {
        i++;
      }
    }
    int j = 0;
    for (const auto& rhs_key : rhs_keys_) {
      if (S.GetEvaluator(rhs_key, my_keys_[0].second).IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
          
        // FIX ME AND MAKE THIS LESS HACKY AND CONST CORRECT --etc
        CompositeVector drhs = S.GetDerivativeW<CompositeVector>(
            rhs_key, my_keys_[0].second, wrt_key, wrt_tag, rhs_key);
        // std::cout << "Adding source " << rhs_key << " op to operator";
        // drhs.Print(std::cout);

        for (const auto& comp : drhs) {
          if (AmanziMesh::entity_kind(comp) == AmanziMesh::CELL) {
            auto op_cell = global_op->getOp(i);
            op_cell->diag->assign(*drhs.GetComponent(comp, false));
            op_cell->diag->scale(rhs_scalars_[j]);
            i++;
            
          } else {
            AMANZI_ASSERT(0);
          }
        }
      }
      j++;
    }
    
  } else {
    // off diagonal
    throw("OperatorApply not implemented yet for offdiagonals!");
  }

  // assemble
  //global_op->AssembleMatrix();
  global_op->initializeInverse();
}


} // namespace Amanzi
