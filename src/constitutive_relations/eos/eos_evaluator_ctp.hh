/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluatorTP is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_CTP_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_CTP_HH_

#include "eos_sw.hh"
#include "Factory.hh"
#include "eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

class EOSEvaluatorCTP : public EOSEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  EOSEvaluatorCTP(Teuchos::ParameterList& plist);

  EOSEvaluatorCTP(const EOSEvaluatorCTP& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone()  const override;

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key,
                                               const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

  Teuchos::RCP<EOS_SW> get_EOS() { return eos_sw_; }
 protected:
  // the actual model
  Teuchos::RCP<EOS_SW> eos_sw_;


  // Keys for fields
  // dependencies
  Key temp_key_;
  Key pres_key_;
  Key conc_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EOSEvaluatorCTP> factory_;
};

} // namespace
} // namespace

#endif
