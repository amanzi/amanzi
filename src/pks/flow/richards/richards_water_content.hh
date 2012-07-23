/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldModel for water content.

Wrapping this conserved quantity as a field model makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in Richards equation.
----------------------------------------------------------------------------- */


#ifndef AMANZI_RICHARDS_WATER_CONTENT_HH_
#define AMANZI_RICHARDS_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "secondary_variable_field_model.hh"

namespace Amanzi {
namespace Flow {

class RichardsWaterContent : public SecondaryVariableFieldModel {

 public:
  RichardsWaterContent(Teuchos::ParameterList& wc_plist, const Teuchos::Ptr<State>& S);
  RichardsWaterContent(const RichardsWaterContent& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

};

#endif
