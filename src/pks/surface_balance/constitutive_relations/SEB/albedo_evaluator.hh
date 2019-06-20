/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (coonet @ ornl.gov)
   
 ------------------------------------------------------------------------- */

//! AlbedoEvaluator: evaluates albedos and emissivities

/*!
Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for two channels -- water/ice/land and
snow.  Note this internally calculates albedo of snow based upon snow density.

Channels are: 0 = land/ice/water, 1 = snow.

* `"albedo ice [-]`" ``[double]`` **0.44** 
* `"albedo water [-]`" ``[double]`` **0.1168** 
* `"albedo ground surface [-]`" ``[double]`` **0.135** Defaults to that of tundra.

* `"emissivity ice [-]`" ``[double]`` **0.98** 
* `"emissivity water [-]`" ``[double]`` **0.995** 
* `"emissivity ground surface [-]`" ``[double]`` **0.92** Defaults to that of tundra.
* `"emissivity snow [-]`" ``[double]`` **0.98**

* `"snow density key`" ``[string]`` **SNOW_DOMAIN-density** 
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth** 
* `"unfrozen fraction key`" ``[string]`` **DOMAIN-unfrozen_fraction**

*/


#ifndef ALBEDO_EVALUATOR_HH_
#define ALBEDO_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

class AlbedoEvaluator : public SecondaryVariablesFieldEvaluator {
 public:
  explicit
  AlbedoEvaluator(Teuchos::ParameterList& plist);
  AlbedoEvaluator(const AlbedoEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new AlbedoEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Key domain_;
  Key domain_snow_;

  Key albedo_key_, emissivity_key_;
  Key snow_dens_key_, ponded_depth_key_, unfrozen_fraction_key_;

  double a_ice_, a_water_, a_tundra_;
  double e_snow_, e_ice_, e_water_, e_tundra_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,AlbedoEvaluator> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
