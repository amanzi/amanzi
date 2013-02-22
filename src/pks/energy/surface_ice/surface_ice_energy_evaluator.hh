/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for internal energy.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

IE = phi * (s_liquid * n_liquid * u_liquid + s_gas * n_gas * u_gas )
  + (1 - phi) * rho_rock * u_rock

This is simply the conserved quantity in the energy equation.
----------------------------------------------------------------------------- */


#ifndef AMANZI_SURFACE_ICE_ENERGY_EVALUATOR_HH_
#define AMANZI_SURFACE_ICE_ENERGY_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {

class SurfaceIceEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SurfaceIceEnergyEvaluator(Teuchos::ParameterList& energy_plist);
  SurfaceIceEnergyEvaluator(const SurfaceIceEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key dens_key_;
  Key dens_ice_key_;
  Key ie_key_;
  Key ie_ice_key_;
  Key uf_key_;
  Key height_key_;

};

} // namespace
} // namespace

#endif
