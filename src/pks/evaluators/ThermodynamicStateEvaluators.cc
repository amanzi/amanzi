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
#include "IAPWS97.hh"
#include "MeshFactory.hh"
#include "PK_Physical.hh"
#include "State.hh"
#include "VerboseObject.hh"
#include "ThermodynamicStateEvaluators.hh"

namespace Amanzi {
namespace Evaluators {

// conversion factor from kJ/kg to J/mol
static constexpr double cfactor = 1000.0 * CommonDefs::MOLAR_MASS_H2O;

/* ******************************************************************
* Thermodynamic state
****************************************************************** */
ThermodynamicStateEvaluator::ThermodynamicStateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("thermodynamic_state", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  enthalpy_key_ = prefix + "enthalpy";

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));

  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
ThermodynamicStateEvaluator::ThermodynamicStateEvaluator(const ThermodynamicStateEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_),
    eos_(other.eos_)
{};


Teuchos::RCP<Evaluator>
ThermodynamicStateEvaluator::Clone() const
{
  return Teuchos::rcp(new ThermodynamicStateEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
ThermodynamicStateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& h_c = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double p, T, v, ap, av, bp, cp, cv, kt;
  for (int c = 0; c != ncells; ++c) {
    double pMPa = p_c[0][c] * 1.0e-6;
    double hkJ = h_c[0][c] / cfactor;
 
    AmanziEOS::Properties prop, liquid, vapor;
    try {
      std::tie(prop, liquid, vapor) = eos_->ThermodynamicsPH(pMPa, hkJ);
    } catch (...) {
      AMANZI_ASSERT(0);
    }

    result_v[(int)TS_t::RGN][c] = prop.rgn;
    result_v[(int)TS_t::T][c] = prop.T;
    result_v[(int)TS_t::RHO][c] = prop.rho;
    result_v[(int)TS_t::V][c] = prop.v;
    result_v[(int)TS_t::CP][c] = prop.cp * 1.0e+3;
    result_v[(int)TS_t::CV][c] = prop.cv * 1.0e+3;
    result_v[(int)TS_t::KT][c] = prop.kt * 1.0e-6;
    result_v[(int)TS_t::AV][c] = prop.av;
    result_v[(int)TS_t::AP][c] = prop.ap;
    result_v[(int)TS_t::BP][c] = prop.bp;
    result_v[(int)TS_t::K][c] = (prop.x == 0.0) ? liquid.k : vapor.k;
    result_v[(int)TS_t::MU][c] = (prop.x == 0.0) ? liquid.mu : vapor.mu;

    // vapor extension
    result_v[(int)TS_t::VV][c] = vapor.v;
    result_v[(int)TS_t::X][c] = prop.x;

    p = p_c[0][c];
    T = prop.T;
    v = prop.v;

    if (prop.rgn == 3) { // Helmholtz
      ap = prop.ap;
      bp = prop.bp;
      cv = prop.cv * 1.0e+3;

      result_v[(int)TS_t::dRHOdP][c] = (cv + p * v * ap) / (v * v * p * (bp * cv + p * ap * ap * T));
      result_v[(int)TS_t::dRHOdH][c] = -ap / (v * v * (cv * bp + p * T * ap * ap)) / CommonDefs::MOLAR_MASS_H2O;

      result_v[(int)TS_t::dTdP][c] = (T * ap - v * bp) / (bp * cv + p * T * ap * ap);
      result_v[(int)TS_t::dTdH][c] = bp / (p * T * ap * ap + bp * cv) / CommonDefs::MOLAR_MASS_H2O;
    }
    else if (prop.rgn == 4) {
      double dvdp_l, dvdp_v, dhdp_l, dhdp_v; // full derivatives
      double x, vl, vv, hl, hv, tmp1, tmp2, dvdh;

      vl = liquid.v;
      hl = liquid.h * 1.0e+3;

      vv = vapor.v;
      hv = vapor.h * 1.0e+3;

      dvdh = (vv - vl) / (hv - hl);
      
      av = liquid.av;
      cp = liquid.cp * 1.0e+3;
      kt = liquid.kt * 1.0e-6;

      dvdp_l = -vl * kt + v * av * T * dvdh;
      dhdp_l = vl * (1 - av * T) + cp * T * dvdh;

      av = vapor.av;
      cp = vapor.cp * 1.0e+3;
      kt = vapor.kt * 1.0e-6;

      dvdp_v = -vv * kt + v * av * T * dvdh;
      dhdp_v = vv * (1 - av * T) + cp * T * dvdh;

      // mechanical (tmp1) and phase change (tmp2) parts
      // this could be re-grouped and combied with identity to cancel h' terms
      x = prop.x;
      tmp1 = (1 - x) * dvdp_l + x * dvdp_v;
      tmp2 = ((1 - x) * dhdp_l + x * dhdp_v) * dvdh; 
      result_v[(int)TS_t::dRHOdP][c] = -(tmp1 - tmp2) / v / v;
      result_v[(int)TS_t::dRHOdH][c] = -dvdh / v / v / CommonDefs::MOLAR_MASS_H2O;

      // Clapeyron
      result_v[(int)TS_t::dTdP][c] = T * dvdh;
      result_v[(int)TS_t::dTdH][c] = 0.0;
 
      result_v[(int)TS_t::K][c] = (1.0 - x) * liquid.k + x * vapor.k;
      result_v[(int)TS_t::MU][c] = (1.0 - x) * liquid.mu + x * vapor.mu;
    }
    else { // Gibbs
      av = prop.av;
      cp = prop.cp * 1.0e+3;
      kt = prop.kt * 1.0e-6;

      result_v[(int)TS_t::dRHOdP][c] = kt / v + av * (1.0 - T * av) / cp;
      result_v[(int)TS_t::dRHOdH][c] = -av / v / cp / CommonDefs::MOLAR_MASS_H2O;

      result_v[(int)TS_t::dTdP][c] = -v * (1.0 - T * av) / cp;
      result_v[(int)TS_t::dTdH][c] = 1.0 / cp / CommonDefs::MOLAR_MASS_H2O;
    }
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void 
ThermodynamicStateEvaluator::EvaluatePartialDerivative_(const State& S,
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
DensityEvaluator::DensityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() != 2) {
    my_keys_.clear();
    my_keys_.push_back(std::make_pair("mass_density_liquid", Tags::DEFAULT));
    my_keys_.push_back(std::make_pair("molar_density_liquid", Tags::DEFAULT));
  }

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  enthalpy_key_ = prefix + "enthalpy";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
DensityEvaluator::DensityEvaluator(const DensityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_)
{}


Teuchos::RCP<Evaluator>
DensityEvaluator::Clone() const
{
  return Teuchos::rcp(new DensityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
DensityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result0_v = *results[0]->ViewComponent("cell");
  auto& result1_v = *results[1]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result0_v[0][c] = ts_c[(int)TS_t::RHO][c];
    result1_v[0][c] = result0_v[0][c] / CommonDefs::MOLAR_MASS_H2O;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
DensityEvaluator::EvaluatePartialDerivative_(const State& S,
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
      result0_v[0][c] = ts_c[(int)TS_t::dRHOdP][c];
      result1_v[0][c] = result0_v[0][c] / CommonDefs::MOLAR_MASS_H2O;
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      result0_v[0][c] = ts_c[(int)TS_t::dRHOdH][c];
      result1_v[0][c] = result0_v[0][c] / CommonDefs::MOLAR_MASS_H2O;
    }
  }
}


/* ******************************************************************
* Temperature evaluator
****************************************************************** */
TemperatureEvaluator::TemperatureEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("temperature", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  enthalpy_key_ = prefix + "enthalpy";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
TemperatureEvaluator::TemperatureEvaluator(const TemperatureEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_)
{}


Teuchos::RCP<Evaluator>
TemperatureEvaluator::Clone() const
{
  return Teuchos::rcp(new TemperatureEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
TemperatureEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = ts_c[(int)TS_t::T][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
TemperatureEvaluator::EvaluatePartialDerivative_(const State& S,
                                                 const Key& wrt_key,
                                                 const Tag& wrt_tag,
                                                 const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double v, p, T, ap, av, bp, cp, cv;
  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ts_c[(int)TS_t::dTdP][c];
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ts_c[(int)TS_t::dTdH][c];
    }
  }
}


/* ******************************************************************
* Thermal conductivity evaluator
****************************************************************** */
ThermalConductivityEvaluator::ThermalConductivityEvaluator(Teuchos::ParameterList& plist)
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

  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
ThermalConductivityEvaluator::ThermalConductivityEvaluator(
    const ThermalConductivityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    density_key_(other.density_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
ThermalConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new ThermalConductivityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
ThermalConductivityEvaluator::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = ts_c[(int)TS_t::K][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
ThermalConductivityEvaluator::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drho, T, dT, tc1, tc2;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  AmanziEOS::Properties prop;

  if (wrt_key == density_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = ts_c[(int)TS_t::T][c];

      if (ts_c[(int)TS_t::RGN][c] == 4.0) {
        double vl, vv, v, x, tcl, tcv;
        v = ts_c[(int)TS_t::V][c];
        x = ts_c[(int)TS_t::X][c];

        vv = ts_c[(int)TS_t::VV][c];
        vl = (v - x * vv) / (1.0 - x);

        tcl = eos_->ThermalConductivity(1.0 / vl, T, prop);
        tcv = eos_->ThermalConductivity(1.0 / vv, T, prop);
        result_v[0][c] = -(tcv - tcl) / (vv - vl) * v * v;
      } else {
        tc1 = ts_c[(int)TS_t::K][c];
        rho = ts_c[(int)TS_t::RHO][c];

        drho = eps * rho;
        tc2 = eos_->ThermalConductivity(rho + drho, T, prop);
        result_v[0][c] = (tc2 - tc1) / drho;
      }
    }

  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = ts_c[(int)TS_t::T][c];
      dT = eps * T;

      if (ts_c[(int)TS_t::RGN][c] == 4.0) {
        double vl, vv, v, x, tcl1, tcl2, tcv1, tcv2;
        v = ts_c[(int)TS_t::V][c];
        x = ts_c[(int)TS_t::X][c];

        vv = ts_c[(int)TS_t::VV][c];
        vl = (v - x * vv) / (1.0 - x);

        tcl1 = eos_->ThermalConductivity(1.0 / vl, T, prop);
        tcl2 = eos_->ThermalConductivity(1.0 / vl, T + dT, prop);

        tcv1 = eos_->ThermalConductivity(1.0 / vv, T, prop);
        tcv2 = eos_->ThermalConductivity(1.0 / vv, T + dT, prop);

        result_v[0][c] = ((1.0 - x) * (tcl2 - tcl1) + x * (tcv2 - tcv1)) / dT;
      } else {
        tc1 = ts_c[(int)TS_t::K][c];
        rho = ts_c[(int)TS_t::RHO][c];
        tc2 = eos_->ThermalConductivity(rho, T + dT, prop);
        result_v[0][c] = (tc2 - tc1) / dT;
      }
    }
  }
}


/* ******************************************************************
* Internal energy liquid evaluator
****************************************************************** */
InternalEnergyLiquidEvaluator::InternalEnergyLiquidEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("temperature", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  enthalpy_key_ = prefix + "enthalpy";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
InternalEnergyLiquidEvaluator::InternalEnergyLiquidEvaluator(
    const InternalEnergyLiquidEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_)
{}


Teuchos::RCP<Evaluator>
InternalEnergyLiquidEvaluator::Clone() const
{
  return Teuchos::rcp(new InternalEnergyLiquidEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
InternalEnergyLiquidEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& h_c = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = h_c[0][c] - p_c[0][c] * ts_c[(int)TS_t::V][c] * CommonDefs::MOLAR_MASS_H2O;
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
InternalEnergyLiquidEvaluator::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& h_c = *S.Get<CompositeVector>(enthalpy_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drhodp, drhodh;
  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      rho = ts_c[(int)TS_t::RHO][c];
      drhodp = ts_c[(int)TS_t::dRHOdP][c];
      result_v[0][c] = -1.0 / rho + drhodp / rho / rho;
    }
  } else if (wrt_key == enthalpy_key_) {
    for (int c = 0; c != ncells; ++c) {
      rho = ts_c[(int)TS_t::RHO][c];
      drhodh = ts_c[(int)TS_t::dRHOdH][c];
      result_v[0][c] = 1.0 + drhodh / rho / rho;
    }
  }
}


/* ******************************************************************
* Isothermal compressibility evaluator 1/rho (drho/dp)_T 
****************************************************************** */
IsothermalCompressibilityEvaluator::IsothermalCompressibilityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair("isothermal_compressibility", Tags::DEFAULT));

  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);
  pressure_key_ = prefix + "pressure";
  enthalpy_key_ = prefix + "enthalpy";

  dependencies_.insert(std::make_pair("thermodynamic_state", Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(enthalpy_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
IsothermalCompressibilityEvaluator::IsothermalCompressibilityEvaluator(
    const IsothermalCompressibilityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    pressure_key_(other.pressure_key_),
    enthalpy_key_(other.enthalpy_key_)
{}


Teuchos::RCP<Evaluator>
IsothermalCompressibilityEvaluator::Clone() const
{
  return Teuchos::rcp(new IsothermalCompressibilityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
IsothermalCompressibilityEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = ts_c[(int)TS_t::KT][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
IsothermalCompressibilityEvaluator::EvaluatePartialDerivative_(
    const State& S,
    const Key& wrt_key,
    const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(false);
}


/* ******************************************************************
* Water/steam viscosity evaluator
****************************************************************** */
ViscosityEvaluator::ViscosityEvaluator(Teuchos::ParameterList& plist)
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

  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));
}


/* ******************************************************************
* Copy operations.
****************************************************************** */
ViscosityEvaluator::ViscosityEvaluator(const ViscosityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    density_key_(other.density_key_),
    temperature_key_(other.temperature_key_)
{}


Teuchos::RCP<Evaluator>
ViscosityEvaluator::Clone() const
{
  return Teuchos::rcp(new ViscosityEvaluator(*this));
}


/* ******************************************************************
* Field value
****************************************************************** */
void
ViscosityEvaluator::Evaluate_(const State& S,
                              const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = ts_c[(int)TS_t::MU][c];
  }
}


/* ******************************************************************
* Field derivative value
****************************************************************** */
void
ViscosityEvaluator::EvaluatePartialDerivative_(const State& S,
                                               const Key& wrt_key,
                                               const Tag& wrt_tag,
                                               const std::vector<CompositeVector*>& results)
{
  const auto& ts_c = *S.Get<CompositeVector>("thermodynamic_state").ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  double rho, drho, T, dT, mu1, mu2;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  if (wrt_key == density_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = ts_c[(int)TS_t::T][c];

      if (ts_c[(int)TS_t::RGN][c] == 4.0) {
        double vl, vv, v, x, mul, muv;
        v = ts_c[(int)TS_t::V][c];
        x = ts_c[(int)TS_t::X][c];

        vv = ts_c[(int)TS_t::VV][c];
        vl = (v - x * vv) / (1.0 - x);

        mul = eos_->Viscosity(1.0 / vl, T);
        muv = eos_->Viscosity(1.0 / vv, T);
        result_v[0][c] = -(muv - mul) / (vv - vl) * v * v;
      } else {
        mu1 = ts_c[(int)TS_t::MU][c];
        rho = ts_c[(int)TS_t::RHO][c];

        drho = eps * rho;
        mu2 = eos_->Viscosity(rho + drho, T);
        result_v[0][c] = (mu2 - mu1) / drho;
      }
    }

  } else if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      T = ts_c[(int)TS_t::T][c];
      dT = eps * T;

      if (ts_c[(int)TS_t::RGN][c] == 4.0) {
        double vl, vv, v, x, mul1, mul2, muv1, muv2;
        v = ts_c[(int)TS_t::V][c];
        x = ts_c[(int)TS_t::X][c];

        vv = ts_c[(int)TS_t::VV][c];
        vl = (v - x * vv) / (1.0 - x);

        mul1 = eos_->Viscosity(1.0 / vl, T);
        mul2 = eos_->Viscosity(1.0 / vl, T + dT);

        muv1 = eos_->Viscosity(1.0 / vv, T);
        muv2 = eos_->Viscosity(1.0 / vv, T + dT);
        result_v[0][c] = ((1.0 - x) * (mul2 - mul1) + x * (muv2 - muv1)) / dT;
      } else {
        mu1 = ts_c[(int)TS_t::MU][c];
        rho = ts_c[(int)TS_t::RHO][c];
        mu2 = eos_->Viscosity(rho, T + dT);
        result_v[0][c] = (mu2 - mu1) / dT;
      }
    }
  }
}

} // namespace Evaluators
} // namespace Amanzi

