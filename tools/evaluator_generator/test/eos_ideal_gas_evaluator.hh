/*
  The ideal gas equation of state evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    modelInitializeParamsList =   cv_ = plist.get<double>("heat capacity");
    namespaceCaps = GENERAL
    evalNameString = ideal gas equation of state
    myMethodArgs = temp_v[0][i], pres_v[0][i]
    myMethodDeclarationArgs = double temp, double pres
    myKey = density
    myKeyFirst = density
    namespace = General
    evalClassName = EosIdealGas
    paramDeclarationList =   double cv_;
    modelMethodDeclaration =   double Density(double temp, double pres) const;
    evalNameCaps = EOS_IDEAL_GAS
    myKeyMethod = Density
    evalName = eos_ideal_gas
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_GENERAL_EOS_IDEAL_GAS_EVALUATOR_HH_
#define AMANZI_GENERAL_EOS_IDEAL_GAS_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace General {
namespace Relations {

class EosIdealGasModel;

class EosIdealGasEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  EosIdealGasEvaluator(Teuchos::ParameterList& plist);
  EosIdealGasEvaluator(const EosIdealGasEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<EosIdealGasModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Key pres_key_;

  Teuchos::RCP<EosIdealGasModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EosIdealGasEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif