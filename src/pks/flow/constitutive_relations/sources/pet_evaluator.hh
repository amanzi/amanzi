/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Ahmad Jan (jana@ornl.gov)
*/

//! Evaluates potential evapotranpiration (PET) 
//! Models are provide in the PERM-IV, Version 4, see pages 90-93, Equations 1-57 to 1-60
/*!

Requires the following dependencies:

* `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
* `"relative humicity key`" ``[string]`` **DOMAIN-relative_humidity**
* `"elevation key`" ``[string]`` **DOMAIN-elevation**
* `"shortwave radiation key`" ``[string]`` **DOMAIN-shortwave_radiation**

*/

#ifndef AMANZI_FLOW_RELATIONS_PET_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_PET_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class PETEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PETEvaluator(Teuchos::ParameterList& plist);
  PETEvaluator(const PETEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new PETEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Exceptions::amanzi_throw("NotImplemented: PETEvaluator currently does not provide derivatives.");
  }

 protected:

  Key at_key_, at_inter_key_, rel_hum_key_, elev_key_, swr_key_;
  double pt_alpha_;
  double HeatFluxDensity(double temp_inter, double temp_next);
  double VPSlope(double temp);
  double PsychrometricConstant(double lh_vap, double evel) ;
  double LatentHeatVaporization(double temp);
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,PETEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
#endif
