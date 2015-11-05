/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_ice * n_ice + s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in the permafrost-Richards equation.
----------------------------------------------------------------------------- */


#ifndef AMANZI_PERMAFROST_WATER_CONTENT_HH_
#define AMANZI_PERMAFROST_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class PermafrostWaterContent : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PermafrostWaterContent(Teuchos::ParameterList& wc_plist);
  PermafrostWaterContent(const PermafrostWaterContent& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key phi_key_;
  Key sl_key_;
  Key mdl_key_;
  Key si_key_;
  Key mdi_key_;
  Key sg_key_;
  Key mdg_key_;
  Key mfg_key_;
  Key cv_key_;

private:
  static Utils::RegisteredFactory<FieldEvaluator,PermafrostWaterContent> reg_;

};

} // namespace
} // namespace
} // namespace

#endif
