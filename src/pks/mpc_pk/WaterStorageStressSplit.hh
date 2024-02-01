/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  MPC PK

  Field evaluator for water storage with a correction for fixed
  stress split.
*/


#ifndef AMANZI_MPC_WATER_STORAGE_STRESS_SPLIT_HH_
#define AMANZI_MPC_WATER_STORAGE_STRESS_SPLIT_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
#include "Tag.hh"

namespace Amanzi {

// Optimal stabilizty coefficient due to Makilic, Wheeler.
// Convergence of iterative coupling for coupled flow and geomechanics.
// Comput Geosci 2013.
inline double
FixedStressStability(double E, double nu, double b)
{
  double mu = E / (2 * (1 + nu));
  return b * b / mu / 2;
}

class WaterStorageStressSplit
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit WaterStorageStressSplit(Teuchos::ParameterList& plist);
  WaterStorageStressSplit(const WaterStorageStressSplit& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  Key saturation_key_, porosity_key_, mol_density_liquid_key_, pressure_key_;
  Key young_modulus_key_, poisson_ratio_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, WaterStorageStressSplit> reg_;
};

} // namespace Amanzi

#endif
