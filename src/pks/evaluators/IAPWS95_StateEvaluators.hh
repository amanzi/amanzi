/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  Collection of evaluators based on IAPWS95 industrial formulation.
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
#include "IAPWS95.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Evaluators {

int constexpr TS95_t_size = 17;
const std::vector<std::string> TS95_names = {
  "mass_density", "enthalpy", "specific_volume", "isobaric_heat_capacity",
  "isocoric_heat_capacity", "isothermal_compressibility", "isobaric_expansion_coef", 
  "relative_pressure_coef", "isothermal_stress_coef", "thermal_conductivity", "viscosity",
  "drhodp", "drhodT", "dudp", "dudT", "vv", "vapor_quality"
};

enum class TS95_t : int {
  RHO = 0,
  H = 1,
  V = 2,
  CP = 3,
  CV = 4,
  KT = 5,
  AV = 6,
  AP = 7,
  BP = 8,
  K = 9,
  MU = 10,
  dRHOdP = 11, // derivatives
  dRHOdT = 12,
  dUdP = 13,
  dUdT = 14,
  VV = 15, // two-phase
  X = 16
};

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
  Key pressure_key_, temperature_key_;
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_StateEvaluator> reg_;
};


class IAPWS95_DensityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
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
  Key pressure_key_, temperature_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_DensityEvaluator> reg_;
};


class IAPWS95_ThermalConductivityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
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
  Key density_key_, temperature_key_;
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
};


class IAPWS95_InternalEnergyEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
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
  Key pressure_key_, temperature_key_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_InternalEnergyEvaluator> reg_;
};


class IAPWS95_ViscosityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
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
  Key density_key_, temperature_key_;
  Teuchos::RCP<AmanziEOS::IAPWS95> eos_;
  static Utils::RegisteredFactory<Evaluator, IAPWS95_ViscosityEvaluator> reg_;
};

} // namespace Evaluators
} // namespace Amanzi

#endif
