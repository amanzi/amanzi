/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Basic interface of a ViscosityModel.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_VISCOSITY_HH_
#define PK_FLOW_VISCOSITY_HH_

#include "secondary_variable_field_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class ViscosityModel : public SecondaryVariableFieldModel {

 public:
  // constructor format for all derived classes
  ViscosityModel(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S);
  ViscosityModel(const ViscosityModel& other);

  // Virtual methods that form the ViscosityModel
  virtual double Viscosity(double T) = 0;
  virtual double DViscosityDT(double T) = 0;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  // PList
  Teuchos::ParameterList visc_plist_;

  // Keys for fields
  // dependencies
  Key temp_key_;
};

} // namespace
} // namespace
} // namespace

#endif
