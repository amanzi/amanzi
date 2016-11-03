/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for converting the darcy flux to volumetric flux

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/
#ifndef AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_
#define AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_

#include "FieldEvaluator_Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class Volumetric_FluxEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  explicit
  Volumetric_FluxEvaluator(Teuchos::ParameterList& plist);

  Volumetric_FluxEvaluator(const Volumetric_FluxEvaluator& other);

  // custom ensure compatibility as all data is not just on the same components
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  Teuchos::RCP<FieldEvaluator> Clone() const;

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Key flux_key_;
  Key dens_key_;
  Key mesh_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,Volumetric_FluxEvaluator> fac_;

};


}//namespace
}//namespace

#endif
