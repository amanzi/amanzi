/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK

  Collection of evaluators based on IAPWS95 formulation for
  pressure/enthalpy primary variables.
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
#include "errors.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "IAPWS95Factory.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"

#include "IAPWS95_StateEvaluatorsPH.hh"

namespace Amanzi {
namespace Evaluators {

/* ******************************************************************
* Thermodynamic state
****************************************************************** */
IAPWS95_StateEvaluatorPH::IAPWS95_StateEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  domain_name_ = plist.template get<std::string>("domain name", "domain");
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_,"thermodynamic_state"), Tags::DEFAULT));

  pressure_key_ = Keys::getKey(domain_name_, "pressure");
  enthalpy_key_ = Keys::getKey(domain_name_, "enthalpy");

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));

  eos_ = AmanziEOS::CreateIAPWS95(plist);
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_StateEvaluatorPH::IAPWS95_StateEvaluatorPH(const IAPWS95_StateEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_),
    eos_(other.eos_)
{};


Teuchos::RCP<Evaluator>
IAPWS95_StateEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new IAPWS95_StateEvaluatorPH(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_StateEvaluatorPH::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  double p, T, rho, v, ap, av, bp, cp, cv, kt;
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& h_c = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    double pMPa = p_c[0][c] * 1.0e-6;
    double hkJ = h_c[0][c] / CommonDefs::ENTHALPY_FACTOR;

    result_c[(int)TSPH_t::dRHOdH][c] = -rho * av / (1000.0 * cp * CommonDefs::MOLAR_MASS_H2O);
 
    AmanziEOS::Properties prop, liquid, vapor;
    try {
      std::tie(prop, liquid, vapor) = eos_->ThermodynamicsPH(pMPa, hkJ);
    } catch (...) {
      Exceptions::amanzi_throw(Errors::CutTimestep());
    }

    result_c[(int)TSPH_t::RHO][c] = prop.rho;
    result_c[(int)TSPH_t::T][c] = prop.T;
    result_c[(int)TSPH_t::V][c] = prop.v;
    result_c[(int)TSPH_t::CP][c] = prop.cp * 1.0e+3;
    result_c[(int)TSPH_t::CV][c] = prop.cv * 1.0e+3;
    result_c[(int)TSPH_t::KT][c] = prop.kt * 1.0e-6;
    result_c[(int)TSPH_t::AV][c] = prop.av;
    result_c[(int)TSPH_t::AP][c] = prop.ap;
    result_c[(int)TSPH_t::BP][c] = prop.bp;
    result_c[(int)TSPH_t::K][c] = (prop.x == 0.0) ? liquid.k : vapor.k;
    result_c[(int)TSPH_t::MU][c] = prop.mu;

    // vapor extension
    result_c[(int)TSPH_t::VV][c] = vapor.v;
    result_c[(int)TSPH_t::X][c] = prop.x;

    p = p_c[0][c];
    v = prop.v;

    T = prop.T;
    rho = prop.rho;
    ap = prop.ap;
    bp = prop.bp;
    cv = prop.cv * 1.0e+3;
    cp = prop.cp * 1.0e+3;

    result_c[(int)TSPH_t::dRHOdP][c] = 1.0 / (v * v * p * bp);
    result_c[(int)TSPH_t::dRHOdH][c] = -rho * av / (1000.0 * cp * CommonDefs::MOLAR_MASS_H2O);

    result_c[(int)TSPH_t::dTdP][c] = 1000.0 * (av * T - 1.0) / (rho * cp);
    result_c[(int)TSPH_t::dTdH][c] = 1.0 / (cp * CommonDefs::MOLAR_MASS_H2O);
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void 
IAPWS95_StateEvaluatorPH::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = 0.0;
  }
}


/* ******************************************************************
* Mass density evaluator
****************************************************************** */
IAPWS95_DensityEvaluatorPH::IAPWS95_DensityEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{

  domain_name_ = plist.template get<std::string>("domain name", "domain");
  if (my_keys_.size() != 2) {
    my_keys_.clear();
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_, "mass_density_liquid"), Tags::DEFAULT));
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_,"molar_density_liquid"), Tags::DEFAULT));
  }

  state_key_ = Keys::getKey(domain_name_, "thermodynamic_state");
  pressure_key_ = Keys::getKey(domain_name_, "pressure");
  enthalpy_key_ = Keys::getKey(domain_name_, "enthalpy");

  dependencies_.insert(std::make_pair(state_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_DensityEvaluatorPH::IAPWS95_DensityEvaluatorPH(const IAPWS95_DensityEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_DensityEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new IAPWS95_DensityEvaluatorPH(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_DensityEvaluatorPH::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  auto& result0_c = *results[0]->ViewComponent("cell");
  auto& result1_c = *results[1]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result0_c[0][c] = ts_c[(int)TSPH_t::RHO][c];
    result1_c[0][c] = result0_c[0][c] / CommonDefs::MOLAR_MASS_H2O;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_DensityEvaluatorPH::EvaluatePartialDerivative_(const State& S,
                                                       const Key& wrt_key,
                                                       const Tag& wrt_tag,
                                                       const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  auto& result0_c = *results[0]->ViewComponent("cell");
  auto& result1_c = *results[1]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result0_c[0][c] = ts_c[(int)TSPH_t::dRHOdP][c];
      result1_c[0][c] = result0_c[0][c] / CommonDefs::MOLAR_MASS_H2O;
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      result0_c[0][c] = ts_c[(int)TSPH_t::dRHOdH][c];
      result1_c[0][c] = result0_c[0][c] / CommonDefs::MOLAR_MASS_H2O;
    }
  }
}


/* ******************************************************************
* Temperature evaluator
****************************************************************** */
IAPWS95_TemperatureEvaluatorPH::IAPWS95_TemperatureEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  domain_name_ = plist.template get<std::string>("domain name", "domain");
  
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_, "temperature"), Tags::DEFAULT));

  //auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = Keys::getKey(domain_name_, "pressure");
  enthalpy_key_ = Keys::getKey(domain_name_, "enthalpy");
  state_key_ = Keys::getKey(domain_name_,"thermodynamic_state");

  dependencies_.insert(std::make_pair(state_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_TemperatureEvaluatorPH::IAPWS95_TemperatureEvaluatorPH(
  const IAPWS95_TemperatureEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_),
    state_key_(other.state_key_),
    domain_name_(other.domain_name_) 
{}


Teuchos::RCP<Evaluator>
IAPWS95_TemperatureEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new IAPWS95_TemperatureEvaluatorPH(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_TemperatureEvaluatorPH::Evaluate_(const State& S,
                                          const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = ts_c[(int)TSPH_t::T][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_TemperatureEvaluatorPH::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double v, p, T, ap, av, bp, cp, cv;
  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ts_c[(int)TSPH_t::dTdP][c];
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ts_c[(int)TSPH_t::dTdH][c];
    }
  }
}


/* ******************************************************************
* Thermal conductivity evaluator
****************************************************************** */
IAPWS95_ThermalConductivityEvaluatorPH::IAPWS95_ThermalConductivityEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  domain_name_ = plist.template get<std::string>("domain name", "domain");
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_,"thermal_conductivity"), Tags::DEFAULT));

  state_key_ = Keys::getKey(domain_name_, "thermodynamic_state");
  density_key_ = Keys::getKey(domain_name_, "mass_density_liquid");
  temperature_key_ = Keys::getKey(domain_name_, "temperature");

  dependencies_.insert(std::make_pair(state_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(density_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  eos_ = AmanziEOS::CreateIAPWS95(plist);
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_ThermalConductivityEvaluatorPH::IAPWS95_ThermalConductivityEvaluatorPH(
    const IAPWS95_ThermalConductivityEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    density_key_(other.density_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_ThermalConductivityEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new IAPWS95_ThermalConductivityEvaluatorPH(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_ThermalConductivityEvaluatorPH::Evaluate_(const State& S,
                                                  const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = ts_c[(int)TSPH_t::K][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_ThermalConductivityEvaluatorPH::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drho, T, dT, tc1, tc2;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  AmanziEOS::Properties prop;

  if (wrt_key == density_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];

      tc1 = ts_c[(int)TSPH_t::K][c];
      rho = ts_c[(int)TSPH_t::RHO][c];

      drho = eps * rho;
      tc2 = eos_->ThermalConductivity(rho + drho, T, prop);
      result_c[0][c] = (tc2 - tc1) / drho;
    }

  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];
      dT = eps * T;

      tc1 = ts_c[(int)TSPH_t::K][c];
      rho = ts_c[(int)TSPH_t::RHO][c];
      tc2 = eos_->ThermalConductivity(rho, T + dT, prop);
      result_c[0][c] = (tc2 - tc1) / dT;
    }
  }
}


/* ******************************************************************
* Internal energy liquid evaluator
****************************************************************** */
IAPWS95_InternalEnergyEvaluatorPH::IAPWS95_InternalEnergyEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  domain_name_ = plist.template get<std::string>("domain name", "domain");
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_,"internal_energy"), Tags::DEFAULT));

  state_key_ = Keys::getKey(domain_name_, "thermodynamic_state");
  pressure_key_ = Keys::getKey(domain_name_, "pressure");
  enthalpy_key_ = Keys::getKey(domain_name_, "enthalpy");

  dependencies_.insert(std::make_pair(state_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_InternalEnergyEvaluatorPH::IAPWS95_InternalEnergyEvaluatorPH(
    const IAPWS95_InternalEnergyEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_InternalEnergyEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new IAPWS95_InternalEnergyEvaluatorPH(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_InternalEnergyEvaluatorPH::Evaluate_(const State& S,
                                             const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& h_c = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = h_c[0][c] - p_c[0][c] * ts_c[(int)TSPH_t::V][c] * CommonDefs::MOLAR_MASS_H2O;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_InternalEnergyEvaluatorPH::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      double rho = ts_c[(int)TSPH_t::RHO][c];
      double drhodp = ts_c[(int)TSPH_t::dRHOdP][c];
      result_c[0][c] = -(1.0 / rho - p_c[0][c] * drhodp / (rho * rho)) * CommonDefs::MOLAR_MASS_H2O;
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      double rho = ts_c[(int)TSPH_t::RHO][c];
      double drhodh = ts_c[(int)TSPH_t::dRHOdH][c];
      result_c[0][c] = 1.0 + 1000.0 * p_c[0][c] * drhodh / (rho * rho);
    }
  }
}


/* ******************************************************************
* Water/steam viscosity evaluator
****************************************************************** */
IAPWS95_ViscosityEvaluatorPH::IAPWS95_ViscosityEvaluatorPH(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  domain_name_ = plist.template get<std::string>("domain name", "domain");
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair(Keys::getKey(domain_name_,"viscosity_liquid"), Tags::DEFAULT));

  state_key_ = Keys::getKey(domain_name_, "thermodynamic_state");
  density_key_ = Keys::getKey(domain_name_, "mass_density_liquid");
  temperature_key_ = Keys::getKey(domain_name_, "temperature");

  dependencies_.insert(std::make_pair(state_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(density_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  eos_ = AmanziEOS::CreateIAPWS95(plist);
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IAPWS95_ViscosityEvaluatorPH::IAPWS95_ViscosityEvaluatorPH(const IAPWS95_ViscosityEvaluatorPH& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    density_key_(other.density_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
IAPWS95_ViscosityEvaluatorPH::Clone() const
{
  return Teuchos::rcp(new IAPWS95_ViscosityEvaluatorPH(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IAPWS95_ViscosityEvaluatorPH::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = ts_c[(int)TSPH_t::MU][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IAPWS95_ViscosityEvaluatorPH::EvaluatePartialDerivative_(const State& S,
                                                         const Key& wrt_key,
                                                         const Tag& wrt_tag,
                                                         const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>(state_key_).ViewComponent("cell");
  const auto& T_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");

  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drho, T, dT, mu1, mu2;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  if (wrt_key == density_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];

      mu1 = ts_c[(int)TSPH_t::MU][c];
      rho = ts_c[(int)TSPH_t::RHO][c];

      drho = eps * rho;
      mu2 = eos_->Viscosity(rho + drho, T);
      result_c[0][c] = (mu2 - mu1) / drho;
    }

  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = T_c[0][c];
      dT = eps * T;

      mu1 = ts_c[(int)TSPH_t::MU][c];
      rho = ts_c[(int)TSPH_t::RHO][c];
      mu2 = eos_->Viscosity(rho, T + dT);
      result_c[0][c] = (mu2 - mu1) / dT;
    }
  }
}

} // namespace Evaluators
} // namespace Amanzi

