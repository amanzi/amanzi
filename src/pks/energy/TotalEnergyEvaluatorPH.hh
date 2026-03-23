/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  FieldEvaluator for the total internal energy. 

    IE = phi * h_liquid * n_liquid + (1 - phi) * rho_rock * u_rock
*/

#ifndef AMANZI_ENERGY_TOTAL_ENERGY_EVALUATOR_PH_HH_
#define AMANZI_ENERGY_TOTAL_ENERGY_EVALUATOR_PH_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {

class TotalEnergyEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit TotalEnergyEvaluatorPH(Teuchos::ParameterList& ee_plist);
  TotalEnergyEvaluatorPH(const TotalEnergyEvaluatorPH& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  Key enthalpy_key_, ie_rock_key_, mol_density_liquid_key_;
  Key particle_density_key_, porosity_key_;
};

} // namespace Energy
} // namespace Amanzi

#endif
