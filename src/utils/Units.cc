/*
  Utils

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/units/base_units/imperial/foot.hpp>
#include <boost/units/base_units/imperial/inch.hpp>

#include "Units.hh"

extern bool Amanzi::Utils::concentration_mol_liter = true;

namespace boost {
namespace units {

std::string symbol_string(const Amanzi::Utils::concentration_amanzi&) {
  return (Amanzi::Utils::concentration_mol_liter) ? "mol/L" : "mol/m^3";
}

std::string name_string(const Amanzi::Utils::concentration_amanzi&) {
  return "concentration";
}

}  // namespace units
}  // namespace boost


namespace Amanzi {
namespace Utils {

namespace bu = boost::units;

/* ******************************************************************
* Constructor creates known units.
****************************************************************** */
void Units::Init(const std::string& conc_units)
{
  pressure_factor_ = 1.0;
  if (conc_units == "molar") {
    concentration_mol_liter = true;
    tcc_factor_ = conversion_factor(bu::si::amount() / liter(), concentration());
  } else {
    concentration_mol_liter = false;
    tcc_factor_ = 1.0;
  }

  // supported units of time (extendable list)
  time_["y"] = 365.25 * 24.0 * 3600.0 * bu::si::second;
  time_["m"] = 365.25 * 2.0 * 3600.0 * bu::si::second;
  time_["d"] = 24.0 * 3600.0 * bu::si::second;
  time_["h"] = 3600.0 * bu::si::second;
  time_["s"] = 1.0 * bu::si::second;

  // supported units of lenght (extendable list)
  length_["km"] = 1000.0 * bu::si::meters;
  length_["m"] = 1.0 * bu::si::meters;
  length_["ft"] = conversion_factor(bu::imperial::foot_base_unit::unit_type(), bu::si::meter) * bu::si::meter;
  length_["in"] = conversion_factor(bu::imperial::inch_base_unit::unit_type(), bu::si::meter) * bu::si::meter;
  length_["cm"] = 0.01 * bu::si::meters;

  // supported units of concentration (extendable list)
  concentration_["mol/m^3"] = 1.0 * concentration();
  concentration_["molar"] = conversion_factor(bu::si::volume(), liter()) * concentration();
  concentration_["ppm"] = 1.0e-3 * concentration();
  concentration_["ppb"] = 1.0e-6 * concentration();
}


/* ******************************************************************
* Convert seconds to any time unit.
****************************************************************** */
double Units::ConvertTime(double t, const std::string& unit)
{
  std::string u(unit);
  boost::algorithm::to_lower(u);

  double val(t);
  if (u == "y" || u == "yr" || u == "year") {
    val /= (365.25 * 24.0 * 3600.0); 
  } else if (u == "d" || u == "day") {
    val /= (24.0 * 3600.0); 
  } else if (u == "m" || u == "month") {
    val /= (365.25 * 2.0 * 3600.0); 
  } else if (u == "h" || u == "hour") {
    val /= 3600.0; 
  }

  return val;
}


/* ******************************************************************
* Convert any input time to any output time.
****************************************************************** */
double Units::ConvertTime(double val,
                          const std::string& in_unit,
                          const std::string& out_unit,
                          bool& flag)
{ 
  flag = true;
  if (time_.find(in_unit) == time_.end() ||
      time_.find(out_unit) == time_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= time_[in_unit].value() / time_[out_unit].value();
  return tmp;
}


/* ******************************************************************
* Convert any input length to any output length.
****************************************************************** */
double Units::ConvertLength(double val,
                            const std::string& in_unit,
                            const std::string& out_unit,
                            bool& flag)
{ 
  flag = true;
  if (length_.find(in_unit) == length_.end() ||
      length_.find(out_unit) == length_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= length_[in_unit].value() / length_[out_unit].value();
  return tmp;
}


/* ******************************************************************
* Convert any input concentration to any output concentration.
****************************************************************** */
double Units::ConvertConcentration(double val,
                                   const std::string& in_unit,
                                   const std::string& out_unit,
                                   double molar_mass,
                                   bool& flag)
{ 
  flag = true;
  if (concentration_.find(in_unit) == concentration_.end() ||
      concentration_.find(out_unit) == concentration_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= concentration_[in_unit].value() / concentration_[out_unit].value();

  // It is not clear how to deal properly with dimentionless units.
  if (in_unit == "ppm" || in_unit == "ppb") tmp /= molar_mass;
  if (out_unit == "ppm" || out_unit == "ppb") tmp *= molar_mass;

  return tmp;
}

}  // namespace Utils
}  // namespace Amanzi

