/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  FieldEvaluator for the total internal energy. Wrapping this conserved
  quantity as a field evaluator makes it easier to take derivatives,
  keep updated, and the like. The equation for this is simply:

    IE = phi * (s_liquid * n_liquid * u_liquid + s_gas * n_gas * u_gas)
       + (1 - phi) * rho_rock * u_rock
*/

#ifndef AMANZI_ENERGY_TOTAL_ENERGY_EVALUATOR_HH_
#define AMANZI_ENERGY_TOTAL_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {

class TotalEnergyEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit TotalEnergyEvaluator(Teuchos::ParameterList& ee_plist);
  TotalEnergyEvaluator(const TotalEnergyEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  bool vapor_diffusion_, include_potential_;
  double liquid_molar_mass_;

  Key particle_density_key_, porosity_key_, sat_liquid_key_;
  Key ie_rock_key_, ie_liquid_key_, ie_gas_key_;
  Key mol_density_liquid_key_, mol_density_gas_key_;

  bool aperture_;
  Key aperture_key_;
};

} // namespace Energy
} // namespace Amanzi

#endif
