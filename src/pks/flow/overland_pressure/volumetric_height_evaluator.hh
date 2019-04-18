/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded depth from ponded depth and subgrid parameters.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.

*/

#ifndef AMANZI_FLOW_RELATIONS_VOLUMETRIC_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_VOLUMETRIC_HEIGHT_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class VolumetricHeightEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  VolumetricHeightEvaluator(Teuchos::ParameterList& plist);
  VolumetricHeightEvaluator(const VolumetricHeightEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new VolumetricHeightEvaluator(*this));
  }

 protected:
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key pd_key_;
  Key delta_max_key_, delta_ex_key_;

  // TODO: put these functions into a model and share them across evaluators:
  //  - this one
  //  - overland_subgrid_water_content_evaluator
  //  - volumetric_ponded_depth evaluator --etc
  double f_(double delta, double del_max, double del_ex) {
    return delta >= del_max ? delta - del_ex:
        std::pow(delta/del_max, 2) * (2*del_max - 3*del_ex)
        + std::pow(delta/del_max,3) * (2*del_ex - del_max);
  }
  double f_prime_(double delta, double del_max, double del_ex) {
    return delta >= del_max ? 1 :
        2 * delta/del_max * (2*del_max - 3*del_ex) / del_max
        + 3 * std::pow(delta/del_max,2) * (2*del_ex - del_max) / del_max;
  }

 private:
  static Utils::RegisteredFactory<FieldEvaluator,VolumetricHeightEvaluator> reg_;

};

} //namespace
} //namespace

#endif
