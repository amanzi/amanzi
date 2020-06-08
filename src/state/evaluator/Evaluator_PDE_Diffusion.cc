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

#include "CSR.hh"
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
  flux_key_ = plist.get<std::string>("flux key");
  my_keys_.emplace_back(std::make_pair(local_op_key_, my_tag_));
  my_keys_.emplace_back(std::make_pair(rhs_key_, my_tag_));
  my_keys_.emplace_back(std::make_pair(flux_key_, my_tag_));

  // dependencies
  // -- tensor
  tensor_coef_key_ = plist.get<std::string>("tensor coefficient key");
  if (tensor_coef_key_ == "none") tensor_coef_key_ = "";
  if (!tensor_coef_key_.empty()) 
    dependencies_.emplace_back(std::make_pair(tensor_coef_key_, my_tag_));

  // -- scalar  
  scalar_coef_key_ = plist.get<std::string>("scalar coefficient key");
  if (scalar_coef_key_ == "none") scalar_coef_key_ = "";
  if (!scalar_coef_key_.empty()) {
    dependencies_.emplace_back(std::make_pair(scalar_coef_key_, my_tag_));
    if (!plist.isParameter("scalar coefficient type")) {
      Errors::Message msg("Evaluator_PDE_Diffusion: \"scalar coefficient type\" not set in list, but key provided.  Please set \"scalar coefficient type\"");
      throw(msg);
    } else if (plist.get<std::string>("scalar coefficient type") == "none") {
      Errors::Message msg("Evaluator_PDE_Diffusion: \"scalar coefficient type\" set to \"none\" but scalar coefficient key provided.  Set a type or remove the key.");
      throw(msg);
    }
  } else {
    plist.set("scalar coefficient type", "none");
  }

  // -- boundary conditions
  bcs_key_ = plist.get<std::string>("boundary conditions key");
  dependencies_.emplace_back(std::make_pair(bcs_key_, my_tag_));

  // -- the vector to which this operator is applied
  //   this depends upon the primary variable directly only if using NLFV or if
  //   we want the Newton correction.  Note we need this name, even if it isn't
  //   a dependency, to check how we are getting derivatives.
  u_key_ = plist.get<std::string>("operator argument key");
  if (plist.get<std::string>("discretization primary") == "nlfv: default" ||
      plist.get<std::string>("Newton correction", "none") != "none") {
    dependencies_.emplace_back(std::make_pair(u_key_, my_tag_));
  }

  // -- density, gravity
  gravity_ = plist.get<bool>("gravity", false);
  if (gravity_) {
    rho_key_ = plist.get<std::string>("density key");
    dependencies_.emplace_back(std::make_pair(rho_key_, my_tag_));
  }
}

void
Evaluator_PDE_Diffusion::EnsureCompatibility(State& S)
{
  // require gravity
  if (gravity_) S.Require<AmanziGeometry::Point>("gravity", "");
  
  // require the rhs
  auto& rhs_fac = S.Require<CompositeVector, CompositeVectorSpace>(
    rhs_key_, my_tag_, rhs_key_);

  // require flux
  auto& flux_fac = S.Require<CompositeVector, CompositeVectorSpace>(
      flux_key_, my_tag_, flux_key_);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  if (rhs_fac.Mesh().get()) {
    mesh = rhs_fac.Mesh();
    flux_fac.SetMesh(mesh)->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  } else if (flux_fac.Mesh().get()) {
    mesh = flux_fac.Mesh();
    rhs_fac.SetMesh(mesh);
  }
  
  if (mesh.get()) {
    // we have a mesh for the RHS, so we can create a diffusion op to get the
    // schema
    Operators::PDE_DiffusionFactory diff_fac;
    if (S.HasDerivativeSet(local_op_key_, my_tag_) &&
        S.GetDerivativeSet(local_op_key_, my_tag_).HasRecord(Keys::getKeyTag(u_key_, my_tag_))) {
      plist_.set("Newton correction", true);
      newton_correction_primary_ = true;
    }
    if (gravity_) {
      pde_ = diff_fac.CreateWithGravity(plist_, mesh);
    } else {
      pde_ = diff_fac.Create(plist_, mesh);
    }
    
    // now we can set up the local op
    auto& lop_fac = S.Require<Operators::Op, Operators::Op_Factory>(
      local_op_key_, my_tag_, local_op_key_);
    lop_fac.set_mesh(mesh);
    Operators::Schema schema(pde_->schema_dofs());
    lop_fac.set_schema(schema);

    // push schema to the rhs cvs
    CompositeVectorSpace cvs = Operators::cvsFromSchema(schema, mesh, true);
    rhs_fac.Update(cvs);

    // require scalar coef on the space required by scalar_coefficient option of
    // operator
    if (!scalar_coef_key_.empty()) {
      S.Require<CompositeVector, CompositeVectorSpace>(scalar_coef_key_, my_tag_)
          .Update(pde_->scalar_coefficient_space());

      // check for derivatives
      if (S.HasDerivativeSet(local_op_key_, my_tag_)) {
        for (const auto& deriv : S.GetDerivativeSet(local_op_key_, my_tag_)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          bool wrt_is_u = wrt.first == u_key_ && wrt.second == my_tag_;
          
          if (S.RequireEvaluator(scalar_coef_key_, my_tag_)
              .IsDifferentiableWRT(S, wrt.first, wrt.second)) {
            // get a local op to store the preconditioner/jacobian ops
            auto& jac_op_fac = S.RequireDerivative<Operators::Op, Operators::Op_Factory>(
                local_op_key_, my_tag_, wrt.first, wrt.second, local_op_key_);
            jac_op_fac.set_mesh(mesh);
            jac_op_fac.set_schema(Operators::Schema(pde_->schema_jacobian()));
            
            // create a PDE for the derivative if needed
            if (!wrt_is_u) {
              Teuchos::ParameterList tmp(plist_);
              tmp.set("exclude primary terms", true);
              tmp.set("Newton correction", true);
              pde_derivs_[deriv.first] = diff_fac.Create(tmp, mesh);
            }

            // require derivatives of the scalar coefficient
            S.RequireDerivative<CompositeVector, CompositeVectorSpace>(scalar_coef_key_, my_tag_,
                    wrt.first, wrt.second)
                .Update(pde_->scalar_coefficient_derivative_space());
          }
        }
      }

      // now we can ensurecompatibility on the coefficient
      S.RequireEvaluator(scalar_coef_key_, my_tag_).EnsureCompatibility(S);
    }

    // require bcs
    auto& bc_fac =
      S.Require<Operators::BCs, Operators::BCs_Factory>(bcs_key_, my_tag_);
    bc_fac.set_mesh(mesh);
    bc_fac.set_entity_kind(AmanziMesh::FACE);
    bc_fac.set_type(WhetStone::DOF_Type::SCALAR);
    S.RequireEvaluator(bcs_key_, my_tag_).EnsureCompatibility(S);

    // require tensors on cells
    if (!tensor_coef_key_.empty()) {
      auto& K_fac =
          S.Require<TensorVector, TensorVector_Factory>(tensor_coef_key_, my_tag_);
      CompositeVectorSpace K_map;
      K_map.SetMesh(mesh);
      K_map.AddComponent("cell", AmanziMesh::CELL, 1);
      K_fac.set_map(K_map);
      S.RequireEvaluator(tensor_coef_key_, my_tag_).EnsureCompatibility(S);
    }

    // require density on cells
    if (gravity_) {
      auto& rho_fac =
          S.Require<CompositeVector,CompositeVectorSpace>(rho_key_, my_tag_);
      rho_fac.SetMesh(mesh)
          ->AddComponent("cell", AmanziMesh::CELL, 1)
          ->SetGhosted();
      S.RequireEvaluator(rho_key_, my_tag_).EnsureCompatibility(S);
    }
  }
}


bool
Evaluator_PDE_Diffusion::Update(State& S, const Key& request)
{
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Evaluator_PDE_Diffusion " << my_keys_[0].first << " requested by "
              << request << "..." << std::endl;
  }

  std::string my_request = Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second);

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = false;
  if (requests_.empty()) {
    // first time through -- give the PDE its data
    update = true;

    auto A_rhs = S.GetPtrW<CompositeVector>(rhs_key_, my_tag_, rhs_key_);
    auto A_lop = S.GetPtrW<Operators::Op>(local_op_key_, my_tag_, local_op_key_);
    pde_->set_local_op(A_lop);
    pde_->global_operator()->set_rhs(A_rhs);
  }
    
  // Check tensor coefficient
  if (S.GetEvaluator(tensor_coef_key_, my_tag_).Update(S, my_request)) {
    update = true;
    pde_->SetTensorCoefficient(S.GetPtr<TensorVector>(tensor_coef_key_, my_tag_));
  }

  // Check scalar coefficient
  if (S.GetEvaluator(scalar_coef_key_, my_tag_).Update(S, my_request)) {
    update = true;

    // if we're going to get the Jacobian later, it saves work now to set dkdp as well.
    Teuchos::RCP<const CompositeVector> dkdp = Teuchos::null;
    if (newton_correction_primary_) {
      S.GetEvaluator(scalar_coef_key_, my_tag_).UpdateDerivative(S, my_request, u_key_, my_tag_);
      dkdp = S.GetDerivativePtr<CompositeVector>(scalar_coef_key_, my_tag_, u_key_, my_tag_);
    }

    pde_->SetScalarCoefficient(S.GetPtr<CompositeVector>(scalar_coef_key_, my_tag_), dkdp);
  }

  // Check BCs
  if (S.GetEvaluator(bcs_key_, my_tag_).Update(S, my_request)) {
    update = true;
    auto bcs = S.GetPtr<Operators::BCs>(bcs_key_, my_tag_);
    pde_->SetBCs(bcs, bcs);
  }

  // Check gravity, rho
  if (gravity_) {
    if (S.GetEvaluator(rho_key_, my_tag_).Update(S, my_request)) {
      update = true;
      auto rho = S.GetPtr<CompositeVector>(rho_key_, my_tag_);
      pde_->SetDensity(rho);

      auto grav = S.Get<AmanziGeometry::Point>("gravity", "");
      pde_->SetGravity(grav);
    }
  }

  // Check u for NLFV
  if (plist_.get<std::string>("discretization primary") == "nlfv: default") {
    if (S.GetEvaluator(u_key_, my_tag_).Update(S, my_request)) {
      update = true;
    }
  }

  // do the update
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "... updating " << my_keys_[0].first << " value."
                << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    Update_(S);
    requests_.clear();
    requests_.insert(request);
    return true;
  } else {
    // Otherwise, see if we have filled this request already.
    if (requests_.find(request) == requests_.end()) {
      requests_.insert(request);
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << my_keys_[0].first
                  << "... has changed, but no need to update. " << std::endl;
      }
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << my_keys_[0].first << "... has not changed. " << std::endl;
      }
      return false;
    }
  }
}



bool
Evaluator_PDE_Diffusion::UpdateDerivative(State& S,
        const Key& requestor, const Key& wrt_key, const Key& wrt_tag)
{
  AMANZI_ASSERT(IsDependency(S, wrt_key, wrt_tag));
  AMANZI_ASSERT(!scalar_coef_key_.empty()); // may be relaxed, e.g. if rho(u)?
  
  bool wrt_is_u = wrt_key == u_key_ && wrt_tag == my_tag_;

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Algebraic Variable d" << my_keys_[0].first << "_d" << wrt_key
              << " requested by " << requestor << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  Key wrt = Keys::getKeyTag(wrt_key, wrt_tag);
  auto& deriv_request_set = deriv_requests_[wrt];

  bool update = false;
  if (deriv_request_set.empty()) {
    // has not been done once, need to give it the jac op
    update = true;
    auto A_jac_lop = S.GetDerivativePtrW<Operators::Op>(local_op_key_, my_tag_,
            wrt_key, wrt_tag, local_op_key_);
    if (wrt_is_u) {
      pde_->set_jacobian_op(A_jac_lop);
    } else {
      pde_derivs_[wrt]->set_jacobian_op(A_jac_lop);
    }
  }

  // -- must update if our our dependencies have changed, as these affect the
  // partial derivatives
  Key my_request = Key{ "d" } +
                   Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second) +
                   "_d" + Keys::getKeyTag(wrt_key, wrt_tag);

  // I don't believe this should be required.  It may be a lot of extra work to
  //  compute ourselves (this is assembling local matrices)
  //  update |= Update(S, my_request);

  // Check tensors
  if (!tensor_coef_key_.empty()) {
    if (S.GetEvaluator(tensor_coef_key_, my_tag_).Update(S, my_request)) {
      update = true;
      // set the Tensor coef
      if (wrt_is_u) {
        pde_->SetTensorCoefficient(S.GetPtr<TensorVector>(tensor_coef_key_, my_tag_));
      } else {
        pde_derivs_[wrt]->SetTensorCoefficient(S.GetPtr<TensorVector>(tensor_coef_key_, my_tag_));
      }
    }
  }

  // check scalars
  AMANZI_ASSERT(!scalar_coef_key_.empty()); // need to set BCs outside of this if true
  if (!scalar_coef_key_.empty()) {
    auto& scalar_eval = S.GetEvaluator(scalar_coef_key_, my_tag_);
    
    bool update_scalar = scalar_eval.Update(S, my_request);
    bool update_scalar_deriv = scalar_eval.IsDifferentiableWRT(S, wrt_key, wrt_tag) &&
                               scalar_eval.UpdateDerivative(S, my_request, wrt_key, wrt_tag);
    if (update_scalar || update_scalar_deriv) {
      update = true;
      if (wrt_is_u) {
        // set the scalar coef
        pde_->SetScalarCoefficient(S.GetPtr<CompositeVector>(scalar_coef_key_, my_tag_),
                S.GetDerivativePtr<CompositeVector>(scalar_coef_key_, my_tag_, wrt_key, wrt_tag));
      } else {
        pde_derivs_[wrt]->SetScalarCoefficient(S.GetPtr<CompositeVector>(scalar_coef_key_, my_tag_),
                S.GetDerivativePtr<CompositeVector>(scalar_coef_key_, my_tag_, wrt_key, wrt_tag));
      }
    }
  }    

  // Check BCs
  if (S.GetEvaluator(bcs_key_, my_tag_).Update(S, my_request)) {
    update = true;
    auto bcs = S.GetPtr<Operators::BCs>(bcs_key_, my_tag_);
    pde_->SetBCs(bcs, bcs);
  }

  // Check gravity, rho
  if (gravity_) {
    if (S.GetEvaluator(rho_key_, my_tag_).Update(S, my_request)) {
      update = true;
      auto rho = S.GetPtr<CompositeVector>(rho_key_, my_tag_);
      pde_->SetDensity(rho);

      auto grav = S.Get<AmanziGeometry::Point>("gravity", "");
      pde_->SetGravity(grav);
    }
  }

  if (wrt_is_u) {
    // potentially update the main matrices
    update |= Update(S, my_request);
  }

  // Do the update
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  ... updating." << std::endl;
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
        *vo_.os() << "  ... not updating but new to this request." << std::endl;
      }
      deriv_request_set.insert(requestor);
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
  // compute local ops
  pde_->global_operator()->Zero();
  auto u = S.GetPtr<CompositeVector>(u_key_, my_tag_).ptr();
  auto flux = S.GetPtrW<CompositeVector>(flux_key_, my_tag_, flux_key_).ptr();
  pde_->UpdateMatrices(Teuchos::null, u);
  pde_->UpdateFlux(u, flux);
  pde_->ApplyBCs(true, true, true);
}

void
Evaluator_PDE_Diffusion::UpdateDerivative_(State& S, const Key& wrt_key,
                                           const Key& wrt_tag)
{
  // compute local ops
  Teuchos::Ptr<const CompositeVector> u =
    S.GetPtr<CompositeVector>(u_key_, my_tag_).ptr();
  Teuchos::Ptr<const CompositeVector> flux =
    S.GetPtr<CompositeVector>(flux_key_, my_tag_).ptr();

  bool wrt_is_u = wrt_key == u_key_ && wrt_tag == my_tag_;
  if (wrt_is_u) {
    pde_->jacobian_op()->Zero();
    pde_->UpdateMatricesNewtonCorrection(flux, u);
    pde_->ApplyBCsJacobian();
  } else {
    Key wrt = Keys::getKeyTag(wrt_key, wrt_tag);
    pde_derivs_[wrt]->jacobian_op()->Zero();
    pde_derivs_[wrt]->UpdateMatricesNewtonCorrection(flux, u);
    pde_derivs_[wrt]->ApplyBCsJacobian();
  }    
}

} // namespace Amanzi
