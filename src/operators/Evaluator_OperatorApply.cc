/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Evaluator_OperatorApply is a secondary evaluator that calculates r = b - Ax
/*
  State

*/


/*!

This is expected to be the residual equation for a Leaf PK, therefore it lives on a domain.

It is expected to be of the form:

r = rhs_1 + ... rhs_k + rhs_A00 + rhs_A01 + ... rhs_Anm - (A00 x0 + A01 x0 + ... A0n x0 + A11 x1 + ... Anm xn)

Where:

  x0                   | is the primary variable associated with this residual
  x1 ... xn            | are off-diagonal primary variables.
  A0i, rhs_A0i         | are diagonal, local operators and their RHS boundary
                       | conditions/sources.
  Aji, rhs_Aji         | are not-necessarily diagonal, local operators and their
                       | RHSs applied to the xj.
  rhs_1 ... rhs_k      | are arbitrary source terms which CANNOT NOT BE AFFECTED BY
                       | boundary conditions (i.e, for conservation equations
                       | discretized using control volume methods, BCs affect only
                       | faces while sources are on cells.

Note that we can infer some constraints here:

- The domain and range of the A0i must be subsets of the r space.
  Realistically, r's space is the union of the domain and range spaces of the
  A0i and the space of the rhs_k.

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
  op0_keys_ = Keys::readKeys(plist_, domain, "diagonal local operators", &empty);
  {
    Teuchos::Array<Key> defaults(op0_keys_.size());
    int i = 0;
    for (const auto& op_key : op0_keys_) {
      defaults[i] = op_key + "_rhs";
      ++i;
    }
    op0_rhs_keys_ = Keys::readKeys(plist_, domain, "diagonal local operator rhss", &defaults);
  }

  // opi, opi rhs's
  int j = 0;
  for (const auto& x_key : x_keys_) {
    // strip the domain from x_key
    Key x_key_suffix = Keys::getVarName(x_key);
    op_keys_.emplace_back(
      Keys::readKeys(plist_, domain, x_key_suffix + " local operators", &empty));

    Teuchos::Array<Key> defaults(op_keys_[j].size());
    int i = 0;
    for (const auto& op_key : op_keys_[j]) {
      defaults[i] = op_key + "_rhs";
      ++i;
    }
    op_rhs_keys_.emplace_back(
      Keys::readKeys(plist_, domain, x_key_suffix + " local operator rhss", &defaults));
  }

  // other rhss
  rhs_keys_ = Keys::readKeys(plist_, domain, "additional rhss", &empty);
  if (rhs_keys_.size() > 0) {
    rhs_scalars_ = plist_.get<Teuchos::Array<double>>("rhs coefficients");
    if (rhs_keys_.size() != rhs_scalars_.size()) {
      Errors::Message message;
      message << "Evaluator_OperatorApply for " << my_keys_[0].first << ": RHS terms names ("
              << (int)rhs_keys_.size() << ") and scalings (" << (int)rhs_scalars_.size()
              << ") are not the same length";
      throw(message);
    }
  }

  // push dependencies
  dependencies_.insert(std::make_pair(x0_key_, my_keys_[0].second));
  for (const auto& x_key : x_keys_) dependencies_.insert(std::make_pair(x_key, my_keys_[0].second));
  for (const auto& rhs : rhs_keys_) dependencies_.insert(std::make_pair(rhs, my_keys_[0].second));
  for (const auto& op_key : op0_keys_)
    dependencies_.insert(std::make_pair(op_key, my_keys_[0].second));
  for (const auto& op_rhs_key : op0_rhs_keys_)
    dependencies_.insert(std::make_pair(op_rhs_key, my_keys_[0].second));
  for (const auto& op_list : op_keys_)
    for (const auto& op_key : op_list)
      dependencies_.insert(std::make_pair(op_key, my_keys_[0].second));
  for (const auto& rhs_list : op_rhs_keys_)
    for (const auto& rhs_key : rhs_list)
      dependencies_.insert(std::make_pair(rhs_key, my_keys_[0].second));

  // primary entity is the entity on which the operator primarily exists.
  // This is CELL for control volume methods and NODE for Lagrangian FE
  // methods.
  //
  // All of x, residual, and rhs must have this entity.
  primary_entity_ = plist_.get<std::string>("primary entity", "cell");
  primary_entity_kind_ = AmanziMesh::createEntityKind(primary_entity_);
}

void
Evaluator_OperatorApply::EnsureCompatibility(State& S)
{
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_keys_.size() == 1);
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(
    my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);
  for (const auto& dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);

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
    // add the primary
    my_fac.AddComponent(primary_entity_, primary_entity_kind_, 1);
    CompositeVectorSpace my_fac_mesh_only(my_fac);

    // -- set the rhss' mesh, need for primary entity
    for (const auto& rhs_key : op0_rhs_keys_) {
      auto& rhs_fac = S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second);
      rhs_fac.Update(my_fac_mesh_only);
    }

    for (const auto& rhs_list : op_rhs_keys_) {
      for (const auto& rhs_key : rhs_list) {
        auto& rhs_fac =
          S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second);
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
      my_fac.Update(S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second));
    for (const auto& rhs_list : op_rhs_keys_)
      for (const auto& rhs_key : rhs_list)
        my_fac.Update(
          S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second));

    // -- x must be the same as r
    S.Require<CompositeVector, CompositeVectorSpace>(x0_key_, my_keys_[0].second).Update(my_fac);
    S.RequireEvaluator(x0_key_, my_keys_[0].second).EnsureCompatibility(S);

    // NOTE: should we ensure that other x's are a superset of the
    // domain-space of the other ops? Otherwise there is nothing to say about
    // the other x's.  --etc

    // -- other RHSs must be a subset of r's space
    for (const auto& rhs_key : rhs_keys_) {
      S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second)
        .Update(my_fac_mesh_only);
      S.RequireEvaluator(rhs_key, my_keys_[0].second).EnsureCompatibility(S);
      if (!S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second)
             .SubsetOf(my_fac)) {
        Errors::Message message;
        message << "Evaluator_OperatorApply for " << my_keys_[0].first << ": RHS " << rhs_key
                << " space is not a subset of my space.";
        throw(message);
      }
    }

    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(my_keys_[0].first, my_keys_[0].second)) {
        auto wrt = Keys::splitKeyTag(deriv.first.get());
        auto tag = wrt.second;
        AMANZI_ASSERT(tag == my_keys_[0].second);
        AMANZI_ASSERT(wrt.first == x0_key_); // NEED TO IMPLEMENT OFF-DIAGONALS EVENTUALLY --etc
        auto& deriv_fac = S.RequireDerivative<Operators::Operator, Operators::Operator_Factory>(
          my_keys_[0].first, my_keys_[0].second, wrt.first, tag, my_keys_[0].first);
        deriv_fac.set_mesh(my_fac.Mesh());

        // FIX ME TO BE NON_SQUARE FOR OFF DIAGONALS? --etc
        deriv_fac.set_cvs(my_fac);
        auto plist = Teuchos::rcp(new Teuchos::ParameterList());
        deriv_fac.set_plist(Teuchos::rcp(new Teuchos::ParameterList(plist_)));
      }
    }
  }

  // require data for derivatives
  if (has_derivs) {
    for (const auto& deriv : S.GetDerivativeSet(my_keys_[0].first, my_keys_[0].second)) {
      auto wrt = Keys::splitKeyTag(deriv.first.get());
      auto tag = wrt.second;
      AMANZI_ASSERT(tag == my_keys_[0].second);
      AMANZI_ASSERT(wrt.first == x0_key_); // NEED TO IMPLEMENT OFF-DIAGONALS EVENTUALLY --etc
      // quasi-linear operators covered by Update() call --etc
      // jacobian terms are covered by ParameterList jacobian option --etc

      // rhs terms
      for (const auto& rhs_key : rhs_keys_) {
        if (S.GetEvaluator(rhs_key, my_keys_[0].second).IsDifferentiableWRT(S, wrt.first, tag)) {
          // require the derivative, and update its factory to match its value
          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
             rhs_key, my_keys_[0].second, wrt.first, tag)
            .Update(S.Require<CompositeVector, CompositeVectorSpace>(rhs_key, my_keys_[0].second));
          // re-call EnsureCompatibility() to make sure these derivative requirements get processed
          S.RequireEvaluator(rhs_key, my_keys_[0].second).EnsureCompatibility(S);
        }
      }
    }
  }
}

// These do the actual work
void
Evaluator_OperatorApply::Update_(State& S)
{
  auto& result = S.GetW<CompositeVector>(my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);

  const auto& x = S.Get<CompositeVector>(x0_key_, my_keys_[0].second);
  x.ScatterMasterToGhosted();
  result.PutScalarMasterAndGhosted(0.);

  int i = 0;
  for (const auto& op_key : op0_keys_) {
    // create the global operator
    Operators::Operator_Factory global_op_fac;
    global_op_fac.set_mesh(result.Mesh());
    global_op_fac.set_cvs(x.Map(), result.Map());
    auto global_op = global_op_fac.Create();

    // do the apply
    S.Get<Operators::Op>(op_key, my_keys_[0].second).ApplyMatrixFreeOp(&*global_op, x, result);
    // std::cout << "Result after op0 op:" << std::endl;
    // result.Print(std::cout);

    // Ax - b
    result.Update(-1.0, S.Get<CompositeVector>(op0_rhs_keys_[i], my_keys_[0].second), 1.0);
    // std::cout << "Result after op0 rhs:" << std::endl;
    // result.Print(std::cout);
    ++i;
  }

  int j = 0;
  for (const auto& op_list : op_keys_) {
    const auto& xj = S.Get<CompositeVector>(x_keys_[j], my_keys_[0].second);
    xj.ScatterMasterToGhosted();

    int k = 0;
    for (const auto& op_key : op_list) {
      // create the global operator
      Operators::Operator_Factory global_op_fac;
      global_op_fac.set_mesh(xj.Mesh());
      global_op_fac.set_cvs(xj.Map(), result.Map());
      auto global_op = global_op_fac.Create();

      // do the apply
      S.Get<Operators::Op>(op_key, my_keys_[0].second).ApplyMatrixFreeOp(&*global_op, xj, result);
      // std::cout << "Result after offdiagonal op:" << std::endl;
      // result.Print(std::cout);

      // Ax - b
      result.Update(-1.0, S.Get<CompositeVector>(op_rhs_keys_[j][k], my_keys_[0].second), 1.0);
      // std::cout << "Result after offdiagonal op rhs:" << std::endl;
      // result.Print(std::cout);
      ++k;
    }
    ++j;
  }

  // add all the additional RHSs
  j = 0;
  for (const auto& rhs_key : rhs_keys_) {
    result.Update(rhs_scalars_[j], S.Get<CompositeVector>(rhs_key, my_keys_[0].second), 1.0);
    ++j;
    // std::cout << "Result after rhs " << j << ":" << std::endl;
    // result.Print(std::cout);
  }
}


// creates the operator for applying inverses
//
// Things to note:
//  - Update_() has been called, so we have the local operators.
//  - UpdateDerivative(dep, wrt) has been called on all dependencies that depend upon wrt.
//
void
Evaluator_OperatorApply::UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag)
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
        std::cout << "Adding diffusion op to operator" << std::endl;
        global_op->OpPushBack(S.GetPtrW<Operators::Op>(op_key, my_keys_[0].second, op_key));
        if (S.GetEvaluator(op_key, my_keys_[0].second).IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
          global_op->OpPushBack(S.GetDerivativePtrW<Operators::Op>(
            op_key, my_keys_[0].second, wrt_key, wrt_tag, op_key));
        }
      }

      int j = 0;
      for (const auto& rhs_key : rhs_keys_) {
        if (S.GetEvaluator(rhs_key, my_keys_[0].second).IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
          // FIX ME AND MAKE THIS LESS HACKY AND CONST CORRECT --etc
          CompositeVector drhs = S.GetDerivativeW<CompositeVector>(
            rhs_key, my_keys_[0].second, wrt_key, wrt_tag, rhs_key);
          // std::cout << "Adding source op to operator";
          // drhs.Print(std::cout);
          for (const auto& comp : drhs) {
            if (AmanziMesh::createEntityKind(comp) == AmanziMesh::Entity_kind::CELL) {
              auto op_cell = Teuchos::rcp(new Operators::Op_Cell_Cell(rhs_key, drhs.Mesh()));
              // clobber the diag
              *op_cell->diag = *drhs.ViewComponent(comp, false);
              op_cell->diag->Scale(rhs_scalars_[j]);
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
    global_op->SymbolicAssembleMatrix();
    global_op->set_inverse_parameters(plist_);
    global_op->InitializeInverse();
  }

  // assemble
  global_op->AssembleMatrix();
  global_op->ComputeInverse();
}


} // namespace Amanzi
