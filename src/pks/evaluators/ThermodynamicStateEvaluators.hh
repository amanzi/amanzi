/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK

  Collection of evaluators based on IAPWS97 industrial formulation.
 
  NOTE: ThermodynamicStateEvaluator is for the auxiliary field which 
  keeps units of the IAPWS97.
*/

#ifndef AMANZI_EVALUATORS_THERMODYNAMIC_STATE_HH_
#define AMANZI_EVALUATORS_THERMODYNAMIC_STATE_HH_

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

int constexpr TS_t_size = 18;
const std::vector<std::string> TS_names = {
  "region", "temperature", "mass_density", "specific_volume", "isobaric_heat_capacity",
  "isocoric_heat_capacity", "isothermal_compressibility", "isobaric_expansion_coef", 
  "relative_pressure_coef", "isothermal_stress_coef", "thermal_conductivity", "viscosity",
  "drhodp", "drhodh", "dtdp", "dtdh",
  "vv", "vapor_quality"
};

enum class TS_t : int {
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

class ThermodynamicStateEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit ThermodynamicStateEvaluator(Teuchos::ParameterList& plist);
  ThermodynamicStateEvaluator(const ThermodynamicStateEvaluator& other);

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
};


class DensityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit DensityEvaluator(Teuchos::ParameterList& plist);
  DensityEvaluator(const DensityEvaluator& other);

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


class TemperatureEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit TemperatureEvaluator(Teuchos::ParameterList& plist);
  TemperatureEvaluator(const TemperatureEvaluator& other);

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


class ThermalConductivityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit ThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityEvaluator(const ThermalConductivityEvaluator& other);

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


class InternalEnergyLiquidEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit InternalEnergyLiquidEvaluator(Teuchos::ParameterList& plist);
  InternalEnergyLiquidEvaluator(const InternalEnergyLiquidEvaluator& other);

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


class IsothermalCompressibilityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit IsothermalCompressibilityEvaluator(Teuchos::ParameterList& plist);
  IsothermalCompressibilityEvaluator(const IsothermalCompressibilityEvaluator& other);

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


class ViscosityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit ViscosityEvaluator(Teuchos::ParameterList& plist);
  ViscosityEvaluator(const ViscosityEvaluator& other);

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

} // namespace Evaluators
} // namespace Amanzi

#endif
