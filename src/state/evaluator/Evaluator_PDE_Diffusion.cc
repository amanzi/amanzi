/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Wraps a PDE_Diffusion to be an Evaluator.

/*!

Lots of options here, document me!

*/

#include "CSR_Factory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "Op_Factory.hh"
#include "Operator_Factory.hh"

#include "Evaluator_PDE_Diffusion.hh"

namespace Amanzi {

Evaluator_PDE_Diffusion::Evaluator_PDE_Diffusion(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  my_tag_ = plist.get<std::string>("tag");

  // my keys
  rhs_key_ = plist.get<std::string>("rhs key");
  local_op_key_ = plist.get<std::string>("local operator key");
  my_keys_.emplace_back(std::make_pair(local_op_key_, my_tag_));
  my_keys_.emplace_back(std::make_pair(rhs_key_, my_tag_));

  if (plist.get<std::string>("Newton correction", "none") != "none") {
    jac_op_key_ = plist.get<std::string>("local operator Jacobian key");
    my_keys_.emplace_back(std::make_pair(jac_op_key_, my_tag_));
  }

  // dependencies
  tensor_coef_key_ = plist.get<std::string>("tensor coefficient key");
  dependencies_.emplace_back(std::make_pair(tensor_coef_key_, my_tag_));

  scalar_coef_key_ = plist.get<std::string>("scalar coefficient key");
  dependencies_.emplace_back(std::make_pair(scalar_coef_key_, my_tag_));

  bcs_key_ = plist.get<std::string>("boundary conditions key");
  dependencies_.emplace_back(std::make_pair(bcs_key_, my_tag_));

  if (plist_.get<std::string>("discretization primary") == "nlfv: default" ||
      plist.get<std::string>("Newton correction", "none") != "none") {
    u_key_ = plist.get<std::string>("operator argument key");
    dependencies_.emplace_back(std::make_pair(u_key_, my_tag_));
  }
}

void
Evaluator_PDE_Diffusion::EnsureCompatibility(State& S)
{
  // require the rhs
  auto& rhs_fac = S.Require<CompositeVector, CompositeVectorSpace>(
    rhs_key_, my_tag_, rhs_key_);
  if (rhs_fac.Mesh().get()) {
    // we have a mesh for the RHS, so we can create a diffusion op to get the
    // schema
    Operators::PDE_DiffusionFactory diff_fac;
    auto diff = diff_fac.Create(plist_, rhs_fac.Mesh());

    // now we can set up the local op
    auto& lop_fac = S.Require<Operators::Op, Operators::Op_Factory>(
      local_op_key_, my_tag_, local_op_key_);
    lop_fac.set_mesh(rhs_fac.Mesh());
    Operators::Schema schema(diff->schema_dofs());
    lop_fac.set_schema(schema);

    // push schema to the rhs cvs
    CompositeVectorSpace cvs = Operators::cvsFromSchema(schema, rhs_fac.Mesh(), true);
    rhs_fac.Update(cvs);

    // require scalar coef on the space required by little_k option of
    // operator
    S.Require<CompositeVector, CompositeVectorSpace>(scalar_coef_key_, my_tag_)
      .Update(diff->little_k_space());
    S.RequireEvaluator(scalar_coef_key_, my_tag_).EnsureCompatibility(S);

    // require bcs
    auto& bc_fac =
      S.Require<Operators::BCs, Operators::BCs_Factory>(bcs_key_, my_tag_);
    bc_fac.set_mesh(rhs_fac.Mesh());
    bc_fac.set_entity_kind(AmanziMesh::FACE);
    bc_fac.set_type(WhetStone::DOF_Type::SCALAR);
    S.RequireEvaluator(bcs_key_, my_tag_).EnsureCompatibility(S);

    // require tensors on cells
    auto& K_fac =
      S.Require<TensorVector, TensorVector_Factory>(tensor_coef_key_, my_tag_);
    CompositeVectorSpace K_map;
    K_map.SetMesh(rhs_fac.Mesh());
    K_map.AddComponent("cell", AmanziMesh::CELL, 1);
    K_fac.set_map(K_map);
    S.RequireEvaluator(tensor_coef_key_, my_tag_).EnsureCompatibility(S);

    // check and require jacobian
    // FIX ME, make a RequireDerivative(), move to EnsureCompatibleDerivative()
    // --etc
    if (!jac_op_key_.empty()) {
      auto& jac_op_fac = S.Require<Operators::Op, Operators::Op_Factory>(
        jac_op_key_, my_tag_, jac_op_key_);
      jac_op_fac.set_mesh(rhs_fac.Mesh());
      jac_op_fac.set_schema(diff->schema_jacobian());
    }
  }
}

bool
Evaluator_PDE_Diffusion::UpdateDerivative(State& S, const Key& requestor,
                                          const Key& wrt_key,
                                          const Key& wrt_tag)
{
  AMANZI_ASSERT(IsDependency(S, wrt_key, wrt_tag));
  AMANZI_ASSERT(!jac_op_key_.empty());
  bool wrt_is_u = wrt_key == u_key_ && wrt_tag == my_tag_;

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Algebraic Variable d" << my_keys_[0].first << "_d" << wrt_key
              << " requested by " << requestor << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = false;

  // -- must update if our our dependencies have changed, as these affect the
  // partial derivatives
  Key my_request = Key{ "d" } +
                   Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second) +
                   "_d" + Keys::getKeyTag(wrt_key, wrt_tag);
  update |= Update(S, my_request);

  // -- must update if any of our dependencies' derivatives have changed
  // NOTE: some assumptions about what is and is not differentiable
  // -- abs perm not a function of p
  AMANZI_ASSERT(!S.GetEvaluator(tensor_coef_key_, my_tag_)
                   .IsDifferentiableWRT(S, wrt_key, wrt_tag));
  if (!jac_op_key_.empty() && S.GetEvaluator(scalar_coef_key_, my_tag_)
                                .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
    update |= S.GetEvaluator(scalar_coef_key_, my_tag_)
                .UpdateDerivative(S, my_request, wrt_key, wrt_tag);
  }
  if (!u_key_.empty() && !wrt_is_u &&
      S.GetEvaluator(u_key_, my_tag_)
        .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
    update |= S.GetEvaluator(u_key_, my_tag_)
                .UpdateDerivative(S, my_request, wrt_tag, wrt_tag);
  }

  // Do the update
  auto request = std::make_tuple(wrt_key, wrt_tag, requestor);
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  ... updating." << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateDerivative_(S, wrt_key, wrt_tag);
    deriv_requests_.clear();
    deriv_requests_.insert(request);
    return true;
  } else {
    // Otherwise, simply service the request
    if (deriv_requests_.find(request) == deriv_requests_.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "  ... not updating but new to this request." << std::endl;
      }
      deriv_requests_.insert(request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "  ... has not changed." << std::endl;
      }
      return false;
    }
  }
}

void
Evaluator_PDE_Diffusion::Update_(State& S)
{
  // get pointers to the results of this update
  auto A_rhs = S.GetPtrW<CompositeVector>(rhs_key_, my_tag_, rhs_key_);
  auto A_lop = S.GetPtrW<Operators::Op>(local_op_key_, my_tag_, local_op_key_);

  // create the global operator
  Operators::Operator_Factory global_op_fac;
  global_op_fac.set_mesh(A_rhs->Mesh());
  global_op_fac.set_cvs(A_rhs->getMap(), A_rhs->getMap());
  auto global_op = global_op_fac.Create();

  // set up the global operator
  global_op->set_rhs(A_rhs);

  Operators::PDE_DiffusionFactory diff_fac;
  auto pde = diff_fac.Create(plist_, global_op);
  pde->set_local_op(A_lop);

  Teuchos::RCP<const Operators::BCs> bcs =
    S.GetPtr<Operators::BCs>(bcs_key_, my_tag_);
  pde->SetBCs(bcs, bcs);

  Teuchos::RCP<const TensorVector> K =
      S.GetPtr<TensorVector>(tensor_coef_key_, my_tag_);
  pde->SetTensorCoefficient(K);

  // at least this is const!
  Teuchos::RCP<const CompositeVector> kr =
    S.GetPtr<CompositeVector>(scalar_coef_key_, my_tag_);
  pde->SetScalarCoefficient(kr, Teuchos::null);

  // compute local ops
  global_op->Zero();

  Teuchos::Ptr<const CompositeVector> u;
  if (!u_key_.empty()) u = S.GetPtr<CompositeVector>(u_key_, my_tag_).ptr();
  pde->UpdateMatrices(Teuchos::null, u);
  pde->ApplyBCs(true, true, true);
}

void
Evaluator_PDE_Diffusion::UpdateDerivative_(State& S, const Key& wrt_key,
                                           const Key& wrt_tag)
{
  // need the spaces
  const CompositeVector& A_rhs = S.Get<CompositeVector>(rhs_key_, my_tag_);

  // create the global operator
  Operators::Operator_Factory global_op_fac;
  global_op_fac.set_mesh(A_rhs.Mesh());
  global_op_fac.set_cvs(A_rhs.getMap(), A_rhs.getMap());
  auto global_op = global_op_fac.Create();

  // set up the global operator
  Operators::PDE_DiffusionFactory diff_fac;
  Teuchos::ParameterList tmp(plist_);
  tmp.set("exclude primary terms",
          true); // turn off the primary terms, just the Jacobian
  auto pde = diff_fac.Create(plist_, global_op);
  if (!jac_op_key_.empty()) {
    pde->set_jacobian_op(
      S.GetPtrW<Operators::Op>(jac_op_key_, my_tag_, jac_op_key_));
  }

  // set the bcs
  Teuchos::RCP<const Operators::BCs> bcs =
    S.GetPtr<Operators::BCs>(bcs_key_, my_tag_);
  pde->SetBCs(bcs, bcs);

  // set the tensor coef
  Teuchos::RCP<const TensorVector> K =
      S.GetPtr<TensorVector>(tensor_coef_key_, my_tag_);
  pde->SetTensorCoefficient(K);

  // set the scalar coef and derivatives
  Teuchos::RCP<const CompositeVector> kr =
    S.GetPtr<CompositeVector>(scalar_coef_key_, my_tag_);
  Teuchos::RCP<const CompositeVector> dkr = S.GetDerivativePtr<CompositeVector>(
    scalar_coef_key_, my_tag_, wrt_key, wrt_tag);
  pde->SetScalarCoefficient(kr, dkr);

  // compute local ops
  global_op->Zero();
  Teuchos::Ptr<const CompositeVector> u =
    S.GetPtr<CompositeVector>(u_key_, my_tag_).ptr();
  pde->UpdateMatricesNewtonCorrection(Teuchos::null, u);
  pde->ApplyBCs(true, true, true);
}

} // namespace Amanzi
