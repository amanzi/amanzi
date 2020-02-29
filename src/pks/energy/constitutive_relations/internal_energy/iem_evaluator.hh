/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The IEM Evaluator simply calls the IEM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_EVALUATOR_

#include "Factory.hh"
#include "iem.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {

class IEMEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  IEMEvaluator(Teuchos::ParameterList& plist);
  IEMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM>& iem);
  IEMEvaluator(const IEMEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& results);

  Teuchos::RCP<IEM> get_IEM() { return iem_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Teuchos::RCP<IEM> iem_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IEMEvaluator> factory_;

};

} //namespace
} //namespace

#endif
