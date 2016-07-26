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
#include <boost/units/base_units/imperial/pound.hpp>
#include <boost/units/base_units/imperial/yard.hpp>
#include <boost/units/base_units/imperial/gallon.hpp>

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
void Units::Init()
{
  pressure_factor_ = 1.0;
  if (system_.concentration == "molar") {
    concentration_mol_liter = true;
    tcc_factor_ = conversion_factor(bu::si::amount() / liter(), concentration());
  } else {
    concentration_mol_liter = false;
    tcc_factor_ = 1.0;
  }

  // supported units of time (extendable list)
  time_["y"] = 365.25 * 24.0 * 3600.0 * bu::si::second;
  time_["d"] = 24.0 * 3600.0 * bu::si::second;
  time_["h"] = 3600.0 * bu::si::second;
  time_["min"] = 60.0 * bu::si::second;
  time_["s"] = 1.0 * bu::si::second;

  // supported units of mass (extendable list)
  mass_["kg"] = 1.0 * bu::si::kilogram;
  mass_["g"] = 0.001 * bu::si::kilogram;
  mass_["lb"] = conversion_factor(bu::imperial::pound_base_unit::unit_type(), bu::si::kilogram) * bu::si::kilogram;

  // supported units of lenght (extendable list)
  length_["km"] = 1000.0 * bu::si::meters;
  length_["m"] = 1.0 * bu::si::meters;
  length_["yd"] = conversion_factor(bu::imperial::yard_base_unit::unit_type(), bu::si::meter) * bu::si::meter;
  length_["ft"] = conversion_factor(bu::imperial::foot_base_unit::unit_type(), bu::si::meter) * bu::si::meter;
  length_["in"] = conversion_factor(bu::imperial::inch_base_unit::unit_type(), bu::si::meter) * bu::si::meter;
  length_["cm"] = 0.01 * bu::si::meters;

  // supported units of volume (extendable list)
  volume_["m3"] = 1.0 * bu::si::volume();
  volume_["gal"] = conversion_factor(bu::imperial::gallon_base_unit::unit_type(), bu::si::volume()) * bu::si::volume();
  volume_["L"] = conversion_factor(liter(), bu::si::volume()) * bu::si::volume();

  // supported units of concentration (extendable list)
  concentration_["mol/m^3"] = 1.0 * concentration();
  concentration_["molar"] = conversion_factor(bu::si::volume(), liter()) * concentration();
  concentration_["ppm"] = 1.0e-3 * concentration();
  concentration_["ppb"] = 1.0e-6 * concentration();

  // supported derived units (simple map is suffient)
  AtomicUnitForm form("kg", 1, "m", -1, "s", -2);
  derived_["Pa"] = form;
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
* Convert any input mass to any output mass.
****************************************************************** */
double Units::ConvertMass(double val,
                          const std::string& in_unit,
                          const std::string& out_unit,
                          bool& flag)
{ 
  flag = true;
  if (mass_.find(in_unit) == mass_.end() ||
      mass_.find(out_unit) == mass_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= mass_[in_unit].value() / mass_[out_unit].value();
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
                                   double mol_mass,
                                   bool& flag)
{ 
  flag = true;
  if (concentration_.find(in_unit) == concentration_.end() ||
      concentration_.find(out_unit) == concentration_.end()) {
    flag = false;
    return val;
  }

  if (mol_mass <= 0.0) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= concentration_[in_unit].value() / concentration_[out_unit].value();

  // It is not clear how to deal properly with dimentionless units.
  if (in_unit == "ppm" || in_unit == "ppb") tmp /= mol_mass;
  if (out_unit == "ppm" || out_unit == "ppb") tmp *= mol_mass;

  return tmp;
}


/* ******************************************************************
* Convert any derived input unit to compatible output unit.
* Special case, out_unit="SI", leads to simple conversion.
****************************************************************** */
double Units::ConvertUnitD(double val,
                           const std::string& in_unit,
                           const std::string& out_unit,
                           double mol_mass,
                           bool& flag)
{
  // replace known complex/derived units
  AtomicUnitForm aut = ComputeAtomicUnitForm_(in_unit, &flag);
 
  for (std::map<std::string, AtomicUnitForm>::iterator it = derived_.begin(); it != derived_.end(); ++it) {
    aut.replace(it->first, it->second); 
  }

  int ntime(0), nmass(0), nlength(0);
  double tmp(val);
  const UnitData& in_data = aut.data();

  UnitData::const_iterator it;
  for (it = in_data.begin(); it != in_data.end(); ++it) {
    if (time_.find(it->first) != time_.end()) {
      tmp *= std::pow(time_[it->first].value(), it->second);
      ntime += it->second;
    }
    else if (mass_.find(it->first) != mass_.end()) {
      tmp *= std::pow(mass_[it->first].value(), it->second);
      nmass += it->second;
    }
    else if (length_.find(it->first) != length_.end()) {
      tmp *= std::pow(length_[it->first].value(), it->second);
      nlength += it->second;
    } 
    else if (volume_.find(it->first) != volume_.end()) {
      tmp *= std::pow(volume_[it->first].value(), it->second);
      nlength += 3 * it->second;
    } else {
      flag = false;
      return val;
    }
  }

  // ad-hoc fix
  if (out_unit == "SI") return tmp;

  // convert from default (SI) units
  aut = ComputeAtomicUnitForm_(out_unit, &flag);
  const UnitData& out_data = aut.data();

  for (it = out_data.begin(); it != out_data.end(); ++it) {
    if (time_.find(it->first) != time_.end()) {
      tmp /= std::pow(time_[it->first].value(), it->second);
      ntime -= it->second;
    }
    else if (mass_.find(it->first) != mass_.end()) {
      tmp /= std::pow(mass_[it->first].value(), it->second);
      nmass -= it->second;
    } 
    else if (length_.find(it->first) != length_.end()) {
      tmp /= std::pow(length_[it->first].value(), it->second);
      nlength -= it->second;
    } 
    else if (volume_.find(it->first) != volume_.end()) {
      tmp /= std::pow(volume_[it->first].value(), it->second);
      nlength -= 3 * it->second;
    } else {
      flag = false;
      return val;
    }
  }

  // consistency of units
  if (ntime || nmass || nlength) {
    flag = false;
    return val;
  }

  return tmp;
}



/* ******************************************************************
* Convert unit string
****************************************************************** */
std::string Units::ConvertUnitS(const std::string& in_unit,
                                const UnitsSystem& system)
{
  // parse the input string
  bool flag;
  AtomicUnitForm aut = ComputeAtomicUnitForm_(in_unit, &flag);

  for (std::map<std::string, AtomicUnitForm>::iterator it = derived_.begin(); it != derived_.end(); ++it) {
    aut.replace(it->first, it->second); 
  }
  UnitData in_data = aut.data();

  // replace units
  UnitData out_data;
  for (UnitData::iterator it = in_data.begin(); it != in_data.end(); ++it) {
    if (time_.find(it->first) != time_.end()) {
      out_data[system.time] = it->second;
    }
    else if (mass_.find(it->first) != mass_.end()) {
      out_data[system.mass] = it->second;
    }
    else if (length_.find(it->first) != length_.end()) {
      out_data[system.length] = it->second;
    }
    else if (concentration_.find(it->first) != concentration_.end()) {
      out_data[system.concentration] = it->second;
    }
  }

  // create string
  std::string separator("");
  std::stringstream ss;
  for (UnitData::iterator it = out_data.begin(); it != out_data.end(); ++it) {
    ss << separator << it->first;

    int i = it->second;
    if (i != 1) {
      ss << "^" << i;
    }
    separator = "*";
  }

  return ss.str();
}


/* ******************************************************************
* Parse unit
****************************************************************** */
AtomicUnitForm Units::ComputeAtomicUnitForm_(const std::string& unit, bool* flag)
{
  *flag = true;

  const char* copy = unit.c_str();
  char* tmp1 = strcpy(new char[unit.size()], unit.c_str());
  char* tmp2 = strtok(tmp1, "^/*");
  char separator_pref, separator_suff;

  AtomicUnitForm form;
  UnitData& data = form.data(); 
  std::pair<UnitData::iterator, bool> status;

  separator_pref = ' ';
  while (tmp2 != NULL) {
    std::string atomic_unit(tmp2);

    status = data.insert(std::pair<std::string, int>(atomic_unit, 0)); 
    UnitData::iterator it = status.first;
      
    if (separator_pref == ' ') {
      (it->second)++;
    } else if (separator_pref == '*') {
      (it->second)++;
    } else if (separator_pref == '/') {
      (it->second)--;
    } else if (separator_pref == '^') {
      *flag = false;
      break;
    }

    separator_suff = copy[tmp2 - tmp1 + strlen(tmp2)];
    tmp2 = strtok(NULL, "^/*");

    if (separator_suff == '^') {
      int i = std::strtol(tmp2, NULL, 10);

      if (separator_pref == '*') {
        it->second += i - 1;
      } else if (separator_pref == '/') {
        it->second -= i - 1;
      }
      separator_suff = copy[tmp2 - tmp1 + strlen(tmp2)];
      tmp2 = strtok(NULL, "^/*");
    }

    separator_pref = separator_suff;
  }

  delete[] tmp1;
  return form;
}

}  // namespace Utils
}  // namespace Amanzi

