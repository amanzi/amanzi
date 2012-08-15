/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The IEM Evaluator simply calls the IEM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class IEM; // forward declaration

class IEMEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  IEMEvaluator(Teuchos::ParameterList& iem_plist);
  IEMEvaluator(Teuchos::ParameterList& iem_plist, const Teuchos::RCP<IEM>& iem);
  IEMEvaluator(const IEMEvaluator& other);

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& results);

  Teuchos::RCP<IEM> get_IEM() { return iem_; }

 protected:

  Teuchos::ParameterList iem_plist_;
  Teuchos::RCP<IEM> iem_;
};

} //namespace
} //namespace
} //namespace

#endif
