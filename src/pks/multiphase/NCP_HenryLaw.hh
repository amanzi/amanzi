/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Multiphase PK

  Field evaluator for noninear complimentary problem, function G.
*/

#ifndef AMANZI_MULTIPHASE_NCP_HENRY_LAW_HH_
#define AMANZI_MULTIPHASE_NCP_HENRY_LAW_HH_

#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "VaporLiquidFactory.hh"

namespace Amanzi {
namespace Multiphase {

class NCP_HenryLaw : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  NCP_HenryLaw(Teuchos::ParameterList& plist);
  NCP_HenryLaw(const NCP_HenryLaw& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key pressure_gas_key_, mol_density_liquid_key_, temperature_key_;
  std::shared_ptr<AmanziEOS::VaporLiquid> vapor_liquid_;

  static Utils::RegisteredFactory<Evaluator, NCP_HenryLaw> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
