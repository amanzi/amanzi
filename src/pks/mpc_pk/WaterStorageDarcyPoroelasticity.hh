/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MPC PK

  Field evaluator for water storage with a correction for volumetric strain.
*/


#ifndef AMANZI_MPC_WATER_STORAGE_DARCY_POROELASTICITY_HH_
#define AMANZI_MPC_WATER_STORAGE_DARCY_POROELASTICITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
#include "Tag.hh"

namespace Amanzi {

class WaterStorageDarcyPoroelasticity
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit WaterStorageDarcyPoroelasticity(Teuchos::ParameterList& plist);
  WaterStorageDarcyPoroelasticity(const WaterStorageDarcyPoroelasticity& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  Key pressure_key_, specific_storage_key_;
  Key biot_key_, strain_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, WaterStorageDarcyPoroelasticity> reg_;
};

} // namespace Amanzi

#endif
