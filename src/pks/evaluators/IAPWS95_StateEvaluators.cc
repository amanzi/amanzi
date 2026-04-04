/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK

  Collection of evaluators based on IAPWS95 formulation.
*/

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
#include "CommonDefs.hh"
#include "CompositeVector.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "IAPWS95.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"

#include "IAPWS95_StateEvaluators.hh"

namespace Amanzi {
namespace Evaluators {

/* ******************************************************************
* Thermodynamic state
****************************************************************** */
IAPWS95_StateEvaluator::IAPWS95_StateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("thermodynamic_state", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  temperature_key_ = prefix + "temperature";

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS95(plist));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_StateEvaluator::IAPWS95_StateEvaluator(const IAPWS95_StateEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    temperature_key_(other.temperature_key_),
    eos_(other.eos_)
{};


Teuchos::RCP<Evaluator>
IAPWS95_StateEvaluator::Clone() const
{
  return Teuchos::rcp(new IAPWS95_StateEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_StateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double p, v, ap, av, bp, cp, cv, kt;
  for (int c = 0; c != ncells; ++c) {
    double pMPa = p_c[0][c] * 1.0e-6;
    double T = T_c[0][c];
 
    AmanziEOS::Properties prop, liquid, vapor;
    try {
      std::tie(prop, liquid, vapor) = eos_->ThermodynamicsPT(pMPa, T);
    } catch (...) {
      AMANZI_ASSERT(0);
    }

    result_v[(int)TS95_t::RHO][c] = prop.rho;
    result_v[(int)TS95_t::H][c] = prop.h;
    result_v[(int)TS95_t::V][c] = prop.v;
    result_v[(int)TS95_t::CP][c] = prop.cp * 1.0e+3;
    result_v[(int)TS95_t::CV][c] = prop.cv * 1.0e+3;
    result_v[(int)TS95_t::KT][c] = prop.kt * 1.0e-6;
    result_v[(int)TS95_t::AV][c] = prop.av;
    result_v[(int)TS95_t::AP][c] = prop.ap;
    result_v[(int)TS95_t::BP][c] = prop.bp;
    result_v[(int)TS95_t::K][c] = (prop.x == 0.0) ? liquid.k : vapor.k;
    // result_v[(int)TS95_t::MU][c] = (prop.x == 0.0) ? liquid.mu : vapor.mu;
    result_v[(int)TS95_t::MU][c] = prop.mu;

    // vapor extension
    result_v[(int)TS95_t::VV][c] = vapor.v;
    result_v[(int)TS95_t::X][c] = prop.x;

    p = p_c[0][c];
    v = prop.v;

    ap = prop.ap;
    bp = prop.bp;
    cv = prop.cv * 1.0e+3;

    result_v[(int)TS95_t::dRHOdP][c] = 1.0 / (v * v * p * bp);
    result_v[(int)TS95_t::dRHOdT][c] = -ap / (v * v * bp);

    result_v[(int)TS95_t::dUdP][c] = (1.0 - T * ap) / bp;
    result_v[(int)TS95_t::dUdT][c] = cv + p * ap * (T * ap - 1.0) / bp;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void 
IAPWS95_StateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& results)
{
  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = 0.0;
  }
}


/* ******************************************************************
* Mass density evaluator
****************************************************************** */
IAPWS95_DensityEvaluator::IAPWS95_DensityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() != 2) {
    my_keys_.clear();
    my_keys_.push_back(std::make_pair("mass_density_liquid", Tags::DEFAULT));
    my_keys_.push_back(std::make_pair("molar_density_liquid", Tags::DEFAULT));
  }

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  temperature_key_ = prefix + "temperature";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_DensityEvaluator::IAPWS95_DensityEvaluator(const IAPWS95_DensityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_DensityEvaluator::Clone() const
{
  return Teuchos::rcp(new IAPWS95_DensityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_DensityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result0_v = *results[0]->ViewComponent("cell");
  auto& result1_v = *results[1]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result0_v[0][c] = ts_c[(int)TS95_t::RHO][c];
    result1_v[0][c] = result0_v[0][c] / CommonDefs::MOLAR_MASS_H2O;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_DensityEvaluator::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  auto& result0_v = *results[0]->ViewComponent("cell");
  auto& result1_v = *results[1]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double v, p, T, ap, av, bp, cv, cp, kt;
  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result0_v[0][c] = ts_c[(int)TS95_t::dRHOdP][c];
      result1_v[0][c] = result0_v[0][c] / CommonDefs::MOLAR_MASS_H2O;
    }
  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      result0_v[0][c] = ts_c[(int)TS95_t::dRHOdT][c];
      result1_v[0][c] = result0_v[0][c] / CommonDefs::MOLAR_MASS_H2O;
    }
  }
}


/* ******************************************************************
* Thermal conductivity evaluator
****************************************************************** */
IAPWS95_ThermalConductivityEvaluator::IAPWS95_ThermalConductivityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("thermal_conductivity", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  density_key_ = prefix + "mass_density_liquid";
  temperature_key_ = prefix + "temperature";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(density_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS95(plist));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_ThermalConductivityEvaluator::IAPWS95_ThermalConductivityEvaluator(
    const IAPWS95_ThermalConductivityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    density_key_(other.density_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_ThermalConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new IAPWS95_ThermalConductivityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_ThermalConductivityEvaluator::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = ts_c[(int)TS95_t::K][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_ThermalConductivityEvaluator::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drho, T, dT, tc1, tc2;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  AmanziEOS::Properties prop;

  if (wrt_key == density_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];

      tc1 = ts_c[(int)TS95_t::K][c];
      rho = ts_c[(int)TS95_t::RHO][c];

      drho = eps * rho;
      tc2 = eos_->ThermalConductivity(rho + drho, T, prop);
      result_v[0][c] = (tc2 - tc1) / drho;
    }

  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];
      dT = eps * T;

      tc1 = ts_c[(int)TS95_t::K][c];
      rho = ts_c[(int)TS95_t::RHO][c];
      tc2 = eos_->ThermalConductivity(rho, T + dT, prop);
      result_v[0][c] = (tc2 - tc1) / dT;
    }
  }
}


/* ******************************************************************
* Internal energy liquid evaluator
****************************************************************** */
IAPWS95_InternalEnergyEvaluator::IAPWS95_InternalEnergyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("internal_energy", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  temperature_key_ = prefix + "temperature";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_InternalEnergyEvaluator::IAPWS95_InternalEnergyEvaluator(
    const IAPWS95_InternalEnergyEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_InternalEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new IAPWS95_InternalEnergyEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_InternalEnergyEvaluator::Evaluate_(const State& S,
                                           const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = (ts_c[(int)TS95_t::H][c] - p_c[0][c] * ts_c[(int)TS95_t::V][c]) * CommonDefs::MOLAR_MASS_H2O;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_InternalEnergyEvaluator::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ts_c[(int)TS95_t::dUdP][c];
    }
  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ts_c[(int)TS95_t::dUdT][c];
    }
  }
}


/* ******************************************************************
* Water/steam viscosity evaluator
****************************************************************** */
IAPWS95_ViscosityEvaluator::IAPWS95_ViscosityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("viscosity_liquid", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  density_key_ = prefix + "mass_density_liquid";
  temperature_key_ = prefix + "temperature";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(density_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS95(plist));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_ViscosityEvaluator::IAPWS95_ViscosityEvaluator(const IAPWS95_ViscosityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    density_key_(other.density_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_ViscosityEvaluator::Clone() const
{
  return Teuchos::rcp(new IAPWS95_ViscosityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_ViscosityEvaluator::Evaluate_(const State& S,
                                      const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = ts_c[(int)TS95_t::MU][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_ViscosityEvaluator::EvaluatePartialDerivative_(const State& S,
                                                       const Key& wrt_key,
                                                       const Tag& wrt_tag,
                                                       const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drho, T, dT, mu1, mu2;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  if (wrt_key == density_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];

      mu1 = ts_c[(int)TS95_t::MU][c];
      rho = ts_c[(int)TS95_t::RHO][c];

      drho = eps * rho;
      mu2 = eos_->Viscosity(rho + drho, T);
      result_v[0][c] = (mu2 - mu1) / drho;
    }

  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];
      dT = eps * T;

      mu1 = ts_c[(int)TS95_t::MU][c];
      rho = ts_c[(int)TS95_t::RHO][c];
      mu2 = eos_->Viscosity(rho, T + dT);
      result_v[0][c] = (mu2 - mu1) / dT;
    }
  }
}

} // namespace Evaluators
} // namespace Amanzi

