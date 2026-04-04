/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  Collection of evaluators based on IAPWS97 industrial formulation.
 
  NOTE: IAPWS97_StateEvaluator is for the auxiliary field which 
  keeps units of the IAPWS97.
*/

#ifndef AMANZI_EVALUATORS_IAPWS97_STATE_HH_
#define AMANZI_EVALUATORS_IAPWS97_STATE_HH_

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
#include "IAPWS97.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Evaluators {

int constexpr TS97_t_size = 18;
const std::vector<std::string> TS97_names = {
  "region", "temperature", "mass_density", "specific_volume", "isobaric_heat_capacity",
  "isocoric_heat_capacity", "isothermal_compressibility", "isobaric_expansion_coef", 
  "relative_pressure_coef", "isothermal_stress_coef", "thermal_conductivity", "viscosity",
  "drhodp", "drhodh", "dtdp", "dtdh",
  "vv", "vapor_quality"
};

enum class TS97_t : int {
  RGN = 0,
  T = 1,
  RHO = 2,
  V = 3,
  CP = 4,
  CV = 5,
  KT = 6,
  AV = 7,
  AP = 8,
  BP = 9,
  K = 10,
  MU = 11,
  dRHOdP = 12, // derivatives
  dRHOdH = 13,
  dTdP = 14,
  dTdH = 15,
  VV = 16, // two-phase
  X = 17
};

class IAPWS97_StateEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_StateEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_StateEvaluator(const IAPWS97_StateEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key pressure_key_, enthalpy_key_;
  Teuchos::RCP<AmanziEOS::IAPWS97> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS97_StateEvaluator> reg_;
};


class IAPWS97_DensityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_DensityEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_DensityEvaluator(const IAPWS97_DensityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key pressure_key_, enthalpy_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS97_DensityEvaluator> reg_;
};


class IAPWS97_TemperatureEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_TemperatureEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_TemperatureEvaluator(const IAPWS97_TemperatureEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key pressure_key_, enthalpy_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS97_TemperatureEvaluator> reg_;
};


class IAPWS97_ThermalConductivityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_ThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_ThermalConductivityEvaluator(const IAPWS97_ThermalConductivityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key density_key_, temperature_key_;
  Teuchos::RCP<AmanziEOS::IAPWS97> eos_;
};


class IAPWS97_InternalEnergyEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_InternalEnergyEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_InternalEnergyEvaluator(const IAPWS97_InternalEnergyEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key pressure_key_, enthalpy_key_;
};


class IAPWS97_IsothermalCompressibilityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_IsothermalCompressibilityEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_IsothermalCompressibilityEvaluator(const IAPWS97_IsothermalCompressibilityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key pressure_key_, enthalpy_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS97_IsothermalCompressibilityEvaluator> reg_;
};


class IAPWS97_ViscosityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IAPWS97_ViscosityEvaluator(Teuchos::ParameterList& plist);
  IAPWS97_ViscosityEvaluator(const IAPWS97_ViscosityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) final {};

 private:
  Key density_key_, temperature_key_;
  Teuchos::RCP<AmanziEOS::IAPWS97> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS97_ViscosityEvaluator> reg_;
};

} // namespace Evaluators
} // namespace Amanzi

#endif
