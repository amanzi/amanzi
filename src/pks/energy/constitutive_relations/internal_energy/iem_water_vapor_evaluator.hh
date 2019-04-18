/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The IEM Evaluator simply calls the IEM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_WATER_VAPOR_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_WATER_VAPOR_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "iem_water_vapor.hh"

namespace Amanzi {
namespace Energy {

class IEMWaterVaporEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  IEMWaterVaporEvaluator(Teuchos::ParameterList& plist);
  IEMWaterVaporEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<IEMWaterVapor>& iem);
  IEMWaterVaporEvaluator(const IEMWaterVaporEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& results);

  Teuchos::RCP<IEMWaterVapor> get_IEM() { return iem_; }

 protected:
  void InitializeFromPlist_();

  Key temp_key_;
  Key mol_frac_key_;
  Teuchos::RCP<IEMWaterVapor> iem_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IEMWaterVaporEvaluator> factory_;

};

} //namespace
} //namespace

#endif
