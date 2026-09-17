/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  Collection of evaluators based on IAPWS95 formulation.
*/

#ifndef AMANZI_EVALUATORS_IAPWS95_STATE_HH_
#define AMANZI_EVALUATORS_IAPWS95_STATE_HH_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CompositeVector.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "EvaluatorSecondaryMonotypeDetached.hh"
#include "IAPWS_Helper.hh"
#include "IAPWS95.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Evaluators {

class IAPWS95_StateEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_StateEvaluator(Teuchos::ParameterList& plist);
  IAPWS95_StateEvaluator(const IAPWS95_StateEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key domain_name_;
  Key pressure_key_, temperature_key_;
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_StateEvaluator> reg_;
};


class IAPWS95_DensityEvaluator : public EvaluatorSecondaryMonotypeDetached {
 public:
  explicit IAPWS95_DensityEvaluator(Teuchos::ParameterList& plist);
  IAPWS95_DensityEvaluator(const IAPWS95_DensityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key domain_name_;
  Key pressure_key_, temperature_key_, state_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_DensityEvaluator> reg_;
};


class IAPWS95_ThermalConductivityEvaluator : public EvaluatorSecondaryMonotypeDetached {
 public:
  explicit IAPWS95_ThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  IAPWS95_ThermalConductivityEvaluator(const IAPWS95_ThermalConductivityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key domain_name_;
  Key density_key_, temperature_key_, state_key_;
  Teuchos::RCP<AmanziEOS::IAPWS> eos_;
};


class IAPWS95_InternalEnergyEvaluator : public EvaluatorSecondaryMonotypeDetached {
 public:
  explicit IAPWS95_InternalEnergyEvaluator(Teuchos::ParameterList& plist);
  IAPWS95_InternalEnergyEvaluator(const IAPWS95_InternalEnergyEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key domain_name_;
  Key pressure_key_, temperature_key_, state_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_InternalEnergyEvaluator> reg_;
};


class IAPWS95_ViscosityEvaluator : public EvaluatorSecondaryMonotypeDetached {
 public:
  explicit IAPWS95_ViscosityEvaluator(Teuchos::ParameterList& plist);
  IAPWS95_ViscosityEvaluator(const IAPWS95_ViscosityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key domain_name_;
  Key density_key_, temperature_key_, state_key_;
  Teuchos::RCP<AmanziEOS::IAPWS> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_ViscosityEvaluator> reg_;
};

} // namespace Evaluators
} // namespace Amanzi

#endif
