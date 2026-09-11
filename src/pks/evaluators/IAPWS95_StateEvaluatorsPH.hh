/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  Collection of evaluators based on IAPWS95 formulation for
  pressure/enthalpy primary variables.
*/

#ifndef AMANZI_EVALUATORS_IAPWS95_STATE_PH_HH_
#define AMANZI_EVALUATORS_IAPWS95_STATE_PH_HH_

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
#include "IAPWS_Helper.hh"
#include "IAPWS95.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Evaluators {

class IAPWS95_StateEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_StateEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_StateEvaluatorPH(const IAPWS95_StateEvaluatorPH& other);

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
  Key pressure_key_, enthalpy_key_;
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_StateEvaluatorPH> reg_;
};


class IAPWS95_DensityEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_DensityEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_DensityEvaluatorPH(const IAPWS95_DensityEvaluatorPH& other);

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
  Key pressure_key_, enthalpy_key_, state_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_DensityEvaluatorPH> reg_;
};


class IAPWS95_TemperatureEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_TemperatureEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_TemperatureEvaluatorPH(const IAPWS95_TemperatureEvaluatorPH& other);

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
  Key pressure_key_, enthalpy_key_, state_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_TemperatureEvaluatorPH> reg_;
};


class IAPWS95_ThermalConductivityEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_ThermalConductivityEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_ThermalConductivityEvaluatorPH(const IAPWS95_ThermalConductivityEvaluatorPH& other);

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
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
};


class IAPWS95_InternalEnergyEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_InternalEnergyEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_InternalEnergyEvaluatorPH(const IAPWS95_InternalEnergyEvaluatorPH& other);

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
  Key pressure_key_, enthalpy_key_, state_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_InternalEnergyEvaluatorPH> reg_;
};


class IAPWS95_ViscosityEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_ViscosityEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_ViscosityEvaluatorPH(const IAPWS95_ViscosityEvaluatorPH& other);

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
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_ViscosityEvaluatorPH> reg_;
};


class IAPWS95_IsothermalCompressibilityEvaluatorPH
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS95_IsothermalCompressibilityEvaluatorPH(Teuchos::ParameterList& plist);
  IAPWS95_IsothermalCompressibilityEvaluatorPH(const IAPWS95_IsothermalCompressibilityEvaluatorPH& other);

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
  Key pressure_key_, enthalpy_key_, state_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_IsothermalCompressibilityEvaluatorPH> reg_;
};


} // namespace Evaluators
} // namespace Amanzi

#endif
