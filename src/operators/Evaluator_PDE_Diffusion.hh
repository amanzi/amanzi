/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Wraps a PDE_Diffusion to be an Evaluator.

/*!

Lots of options here, document me!  
  
*/

#ifndef AMANZI_STATE_EVALUATOR_PDE_DIFFUSION_HH_
#define AMANZI_STATE_EVALUATOR_PDE_DIFFUSION_HH_

#include "Teuchos_ParameterList.hpp"

#include "Evaluator_Factory.hh"
#include "State.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {

class Evaluator_PDE_Diffusion : public EvaluatorSecondary {
 public:
  Evaluator_PDE_Diffusion(Teuchos::ParameterList& plist);

  Evaluator_PDE_Diffusion(const Evaluator_PDE_Diffusion& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new Evaluator_PDE_Diffusion(*this));
  };

  virtual void EnsureCompatibility(State& S) override;

  // jac-op handled in EnsureCompatibility()
  // virtual void EnsureCompatibleDerivative(State &S, const Key& wrt_key, const Tag& wrt_tag) override {};

  virtual bool
  UpdateDerivative(State& S, const Key& requestor, const Key& wrt_key, const Tag& wrt_tag) override;

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return IsDependency(S, wrt_key, wrt_tag) && !jac_op_key_.empty();
  }

 protected:
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override;

 protected:
  Tag my_tag_;
  Key rhs_key_, local_op_key_, jac_op_key_;
  Key tensor_coef_key_, scalar_coef_key_;
  Key bcs_key_;
  Key u_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_PDE_Diffusion> fac_;
};

} // namespace Amanzi

#endif
