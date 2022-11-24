/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for water storage which is the conserved quantity 
  in the Richards equation.

  Wrapping this conserved quantity water storage (WS) as a field evaluator 
  makes it easier to take derivatives, keep updated, and the like.
  The equation for this is simply:

    WS = a * phi * (s_liquid * n_liquid + X_gas * s_gas * n_gas)

  where X_gas is the molar fraction of water in the gas phase,
  a is the optional fracture aperture (default is a = 1), and
  s_liquid + s_gas = 1.
*/


#ifndef AMANZI_FLOW_WATER_STORAGE_EVALUATOR_HH_
#define AMANZI_FLOW_WATER_STORAGE_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
#include "Tag.hh"

namespace Amanzi {
namespace Flow {

class WaterStorage : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit WaterStorage(Teuchos::ParameterList& plist);
  WaterStorage(const WaterStorage& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  void Init_();

 protected:
  bool water_vapor_, aperture_;
  Key saturation_key_, porosity_key_, mol_density_liquid_key_, aperture_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, WaterStorage> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
