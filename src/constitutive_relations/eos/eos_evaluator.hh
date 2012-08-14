/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_HH_

#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class EOSEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  enum EOSMode { EOS_MODE_MASS, EOS_MODE_MOLAR, EOS_MODE_BOTH };

  // constructor format for all derived classes
  explicit
  EOSEvaluator(Teuchos::ParameterList& eos_plist);

  EOSEvaluator(const EOSEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results);

 protected:

  // PList
  Teuchos::ParameterList eos_plist_;

  // the actual model
  Teuchos::RCP<EOS> eos_;
  EOSMode mode_;

  // Keys for fields
  // dependencies
  Key temp_key_;
  Key pres_key_;
  Key a_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EOSEvaluator> factory_;
};

} // namespace
} // namespace

#endif
