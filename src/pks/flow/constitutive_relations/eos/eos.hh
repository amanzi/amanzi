/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Basics of an EOS class... likely isn't specific enough for any given EOS,
  and will need other base models.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _PK_FLOW_EOS_HH_
#define _PK_FLOW_EOS_HH_

#include "secondary_variable_field_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOS : public SecondaryVariableFieldModel {

 public:
  // constructor format for all derived classes
  EOS(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S);
  EOS(const EOS& other);

  // Virtual methods that form the EOS
  virtual double Density(double T, double p) = 0;
  virtual double DDensityDT(double T, double p) = 0;
  virtual double DDensityDp(double T, double p) = 0;

  virtual double molar_mass() = 0;
  virtual bool is_molar_basis() = 0;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  // PList
  Teuchos::ParameterList eos_plist_;

  // Keys for fields
  // dependencies
  Key temp_key_;
  Key pres_key_;
};

} // namespace
} // namespace
} // namespace

#endif
