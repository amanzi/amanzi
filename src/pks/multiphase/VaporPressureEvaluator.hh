/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTIPHASE_VAPOR_PRESSURE_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_VAPOR_PRESSURE_EVALUATOR_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "EOS_SaturatedVaporPressure.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

// Multiphase
#include "WRMmp.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class VaporPressureEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  VaporPressureEvaluator(Teuchos::ParameterList& plist, Teuchos::RCP<WRMmpPartition> wrm);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  Teuchos::RCP<WRMmpPartition> wrm_;
  Teuchos::RCP<AmanziEOS::EOS_SaturatedVaporPressure> svp_;

  std::string temperature_key_, mol_density_liquid_key_;
  std::string saturation_liquid_key_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
