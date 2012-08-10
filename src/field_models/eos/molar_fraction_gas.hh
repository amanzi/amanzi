/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef PK_FLOWRELATIONS_MOLAR_FRACTION_GAS_
#define PK_FLOWRELATIONS_MOLAR_FRACTION_GAS_

#include "secondary_variable_field_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class MolarFractionGas : public SecondaryVariableFieldModel {

 public:
  explicit
  MolarFractionGas(Teuchos::ParameterList& mfg_plist);

  MolarFractionGas(const MolarFractionGas& other);
  virtual Teuchos::RCP<FieldModel> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  // PList
  Teuchos::ParameterList mfg_plist_;
  Key temp_key_;

  Teuchos::RCP<VaporPressureModel> sat_vapor_model_;

 private:
  static Utils::RegisteredFactory<FieldModel,MolarFractionGas> factory_;

};

} //namespace
} //namespace
} //namespace

#endif
