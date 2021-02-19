/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The dynamic subgrid model evaluator gets the subgrid parameters and evolve polygons.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_EXCLUDED_VOL_EVALUATOR_
#define AMANZI_FLOWRELATIONS_EXCLUDED_VOL_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ExcludedVolumeEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  ExcludedVolumeEvaluator(Teuchos::ParameterList& plist);
  ExcludedVolumeEvaluator(const ExcludedVolumeEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;
  
protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Key delta_init_key_,delta_evolve_key_,sg_entity_key_;

private:

  static Utils::RegisteredFactory<FieldEvaluator,ExcludedVolumeEvaluator> reg_;  

};
  
} //namespace
} //namespace 

#endif
