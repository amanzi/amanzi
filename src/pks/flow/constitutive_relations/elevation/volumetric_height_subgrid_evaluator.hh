/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded and snow depths from ponded depth and snow depth.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.
* `"snow depth key`" ``[string]`` **DOMAIN-snow_depth**
         The true height of the water surface.

*/

#ifndef AMANZI_FLOW_RELATIONS_VOLUMETRIC_HEIGHT_SUBGRID_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_VOLUMETRIC_HEIGHT_SUBGRID_EVALUATOR_

#include "secondary_variables_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {


class VolumetricHeightSubgridEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  VolumetricHeightSubgridEvaluator(Teuchos::ParameterList& plist);
  VolumetricHeightSubgridEvaluator(const VolumetricHeightSubgridEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new VolumetricHeightSubgridEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

 protected:

  Key vol_pd_key_, vol_sd_key_;
  Key pd_key_, sd_key_;
  Key delta_max_key_, delta_ex_key_;

  Key domain_snow_, domain_surf_;

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

  bool compatibility_checked_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,VolumetricHeightSubgridEvaluator> reg_;

};

} //namespace
} //namespace

#endif
