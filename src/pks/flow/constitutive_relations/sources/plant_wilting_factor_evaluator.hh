/*
  The plant wilting factor evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).


    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_PLANT_WILTING_FACTOR_EVALUATOR_HH_
#define AMANZI_FLOW_PLANT_WILTING_FACTOR_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class PlantWiltingFactorModel;

class PlantWiltingFactorEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PlantWiltingFactorEvaluator(Teuchos::ParameterList& plist);
  PlantWiltingFactorEvaluator(const PlantWiltingFactorEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  void InitializeFromPlist_();

  Key pc_key_;

  std::vector<std::pair<std::string,Teuchos::RCP<PlantWiltingFactorModel> > > models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PlantWiltingFactorEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
