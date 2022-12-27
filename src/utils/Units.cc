/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Utils

*/

#include <string>

#include "Units.hh"

namespace Amanzi {
namespace Utils {

/* ******************************************************************
* Constructor creates known units.
****************************************************************** */
void
Units::Init()
{
  // supported units of time (extendable list)
  time_["y"] = 365.25 * 24.0 * 3600.0;
  time_["noleap"] = 365. * 24.0 * 3600.0;
  time_["d"] = 24.0 * 3600.0;
  time_["h"] = 3600.0;
  time_["min"] = 60.0;
  time_["s"] = 1.0;

  // supported units of mass (extendable list)
  mass_["ton"] = 1000.0;
  mass_["kg"] = 1.0;
  mass_["g"] = 0.001;
  mass_["lb"] = 0.45359237;

  // supported units of lenght (extendable list)
  length_["km"] = 1000.0;
  length_["m"] = 1.0;
  length_["yd"] = 0.9144;
  length_["ft"] = 0.3048;
  length_["in"] = 0.0254;
  length_["cm"] = 0.01;

  // supported units of volume (extendable list)
  volume_["m3"] = 1.0;
  volume_["gal"] = 0.0037854117840007;
  volume_["L"] = 0.001;

  // supported units of concentration (extendable list)
  concentration_["SI"] = 1.0;
  concentration_["molar"] = 1000.0;
  concentration_["ppm"] = 1.0e-3;
  concentration_["ppb"] = 1.0e-6;

  // supported units of amount of substance (extendable list)
  amount_["mol"] = 1.0;

  // supported units of temperature (extendable list)
  temperature_["K"] = 1.0;
  temperature_["C"] = 1.0;
  temperature_["F"] = 1.0;

  // supported derived units (simple map is suffient)
  derived_["Pa"] = AtomicUnitForm("kg", 1, "m", -1, "s", -2);
  derived_["J"] = AtomicUnitForm("m", 2, "kg", 1, "s", -2);

  // static convertion factor
  concentration_factor_ = concentration_[system_.concentration] / concentration_["SI"];
}


/* ******************************************************************
* Convert any input time to any output time.
****************************************************************** */
double
Units::ConvertTime(double val, const std::string& in_unit, const std::string& out_unit, bool& flag)
{
  flag = true;
  if (time_.find(in_unit) == time_.end() || time_.find(out_unit) == time_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= time_[in_unit] / time_[out_unit];
  return tmp;
}


/* ******************************************************************
* Convert any input mass to any output mass.
****************************************************************** */
double
Units::ConvertMass(double val, const std::string& in_unit, const std::string& out_unit, bool& flag)
{
  flag = true;
  if (mass_.find(in_unit) == mass_.end() || mass_.find(out_unit) == mass_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= mass_[in_unit] / mass_[out_unit];
  return tmp;
}


/* ******************************************************************
* Convert any input length to any output length.
****************************************************************** */
double
Units::ConvertLength(double val,
                     const std::string& in_unit,
                     const std::string& out_unit,
                     bool& flag)
{
  flag = true;
  if (length_.find(in_unit) == length_.end() || length_.find(out_unit) == length_.end()) {
    flag = false;
    return val;
  }

  double tmp(val);
  tmp *= length_[in_unit] / length_[out_unit];
  return tmp;
}


/* ******************************************************************
* Convert any input concentration to any output concentration.
****************************************************************** */
double
Units::ConvertConcentration(double val,
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
  tmp *= concentration_[in_unit] / concentration_[out_unit];

  // It is not clear how to deal properly with dimentionless units.
  if (in_unit == "ppm" || in_unit == "ppb") tmp /= mol_mass;
  if (out_unit == "ppm" || out_unit == "ppb") tmp *= mol_mass;

  return tmp;
}


/* ******************************************************************
* Convert any derived input unit to compatible output unit.
* Special case, out_unit="SI", leads to conversion to SI units.
****************************************************************** */
double
Units::ConvertUnitD(double val,
                    const std::string& in_unit,
                    const std::string& out_unit,
                    double mol_mass,
                    bool& flag)
{
  // replace known complex/derived units
  AtomicUnitForm auf = ComputeAtomicUnitForm_(in_unit, &flag);

  for (std::map<std::string, AtomicUnitForm>::iterator it = derived_.begin(); it != derived_.end();
       ++it) {
    auf.replace(it->first, it->second);
  }

  // replace known atomic units
  if (CompareAtomicUnitForms_(auf, AtomicUnitForm("mol", 1, "L", -1))) {
    auf = AtomicUnitForm("molar", 1);
  } else if (CompareAtomicUnitForms_(auf, AtomicUnitForm("mol", 1, "m", -3))) {
    auf = AtomicUnitForm("SI", 1);
  }

  int ntime(0), nmass(0), nlength(0), nconcentration(0), namount(0), ntemperature(0);
  double tmp(val);
  const UnitData& in_data = auf.data();

  UnitData::const_iterator it;
  for (it = in_data.begin(); it != in_data.end(); ++it) {
    if (time_.find(it->first) != time_.end()) {
      tmp *= std::pow(time_[it->first], it->second);
      ntime += it->second;
    } else if (mass_.find(it->first) != mass_.end()) {
      tmp *= std::pow(mass_[it->first], it->second);
      nmass += it->second;
    } else if (length_.find(it->first) != length_.end()) {
      tmp *= std::pow(length_[it->first], it->second);
      nlength += it->second;
    } else if (volume_.find(it->first) != volume_.end()) {
      tmp *= std::pow(volume_[it->first], it->second);
      nlength += 3 * it->second;
    } else if (concentration_.find(it->first) != concentration_.end()) {
      tmp *= std::pow(concentration_[it->first], it->second);
      if (it->first == "ppm" || it->first == "ppb") tmp /= mol_mass;
      tmp /= concentration_factor_; // for non-SI unit "molar"
      nconcentration += it->second;
    } else if (amount_.find(it->first) != amount_.end()) {
      namount += it->second;
    } else if (temperature_.find(it->first) != temperature_.end()) {
      if (it->first == "C") tmp += 273.15;
      if (it->first == "F") tmp = (tmp + 459.67) * 5.0 / 9;
      ntemperature += it->second;
    } else {
      flag = false;
      return val;
    }
  }

  // ad-hoc fix
  if (out_unit == "SI") return tmp;

  // convert from default (SI) units
  auf = ComputeAtomicUnitForm_(out_unit, &flag);
  const UnitData& out_data = auf.data();

  for (it = out_data.begin(); it != out_data.end(); ++it) {
    if (time_.find(it->first) != time_.end()) {
      tmp /= std::pow(time_[it->first], it->second);
      ntime -= it->second;
    } else if (mass_.find(it->first) != mass_.end()) {
      tmp /= std::pow(mass_[it->first], it->second);
      nmass -= it->second;
    } else if (length_.find(it->first) != length_.end()) {
      tmp /= std::pow(length_[it->first], it->second);
      nlength -= it->second;
    } else if (volume_.find(it->first) != volume_.end()) {
      tmp /= std::pow(volume_[it->first], it->second);
      nlength -= 3 * it->second;
    } else if (concentration_.find(it->first) != concentration_.end()) {
      tmp /= std::pow(concentration_[it->first], it->second);
      if (it->first == "ppm" || it->first == "ppb") tmp *= mol_mass;
      tmp *= concentration_factor_;
      nconcentration -= it->second;
    } else if (amount_.find(it->first) != amount_.end()) {
      namount -= it->second;
    } else if (temperature_.find(it->first) != temperature_.end()) {
      if (it->first == "C") tmp -= 273.15;
      if (it->first == "F") tmp = tmp * 1.8 - 459.67;
      ntemperature -= it->second;
    } else {
      flag = false;
      return val;
    }
  }

  // consistency of units
  if (ntime || nmass || nlength || namount || nconcentration || ntemperature) {
    flag = false;
    return val;
  }

  return tmp;
}


/* ******************************************************************
* Convert unit string
****************************************************************** */
std::string
Units::ConvertUnitS(const std::string& in_unit, const UnitsSystem& system)
{
  // parse the input string
  bool flag;
  AtomicUnitForm auf = ComputeAtomicUnitForm_(in_unit, &flag);

  for (std::map<std::string, AtomicUnitForm>::iterator it = derived_.begin(); it != derived_.end();
       ++it) {
    auf.replace(it->first, it->second);
  }
  UnitData in_data = auf.data();

  // replace units
  UnitData out_data;
  for (UnitData::iterator it = in_data.begin(); it != in_data.end(); ++it) {
    if (time_.find(it->first) != time_.end()) {
      out_data[system.time] = it->second;
    } else if (mass_.find(it->first) != mass_.end()) {
      out_data[system.mass] = it->second;
    } else if (length_.find(it->first) != length_.end()) {
      out_data[system.length] = it->second;
    } else if (concentration_.find(it->first) != concentration_.end()) {
      out_data[system.concentration] = it->second;
    } else if (amount_.find(it->first) != amount_.end()) {
      out_data[system.amount] = it->second;
    } else if (temperature_.find(it->first) != temperature_.end()) {
      out_data[system.temperature] = it->second;
    }
  }

  // create string
  std::string separator("");
  std::stringstream ss;
  for (UnitData::iterator it = out_data.begin(); it != out_data.end(); ++it) {
    ss << separator << it->first;

    int i = it->second;
    if (i != 1) { ss << "^" << i; }
    separator = "*";
  }

  return ss.str();
}


/* ******************************************************************
* Parse unit
****************************************************************** */
AtomicUnitForm
Units::ComputeAtomicUnitForm_(const std::string& unit, bool* flag)
{
  *flag = true;

  const char* copy = unit.c_str();
  char* tmp1 = strcpy(new char[unit.size() + 1], unit.c_str());
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

      if (separator_pref == '*' || separator_pref == ' ') {
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


/* ******************************************************************
* Do two atomic units desribe the same physical quantity?
****************************************************************** */
bool
Units::CompareAtomicUnitForms_(const AtomicUnitForm& auf1, const AtomicUnitForm& auf2)
{
  const UnitData& data1 = auf1.data();
  const UnitData& data2 = auf2.data();

  int ntime(0), nmass(0), nlength(0), nconcentration(0), namount(0), ntemperature(0);
  UnitData::const_iterator it;

  for (it = data1.begin(); it != data1.end(); ++it) {
    if (time_.find(it->first) != time_.end()) ntime += it->second;
    if (mass_.find(it->first) != mass_.end()) nmass += it->second;
    if (length_.find(it->first) != length_.end()) nlength += it->second;
    if (volume_.find(it->first) != volume_.end()) nlength += 3 * it->second;
    if (concentration_.find(it->first) != concentration_.end()) nconcentration += it->second;
    if (amount_.find(it->first) != amount_.end()) namount += it->second;
    if (temperature_.find(it->first) != temperature_.end()) ntemperature += it->second;
  }

  for (it = data2.begin(); it != data2.end(); ++it) {
    if (time_.find(it->first) != time_.end()) ntime -= it->second;
    if (mass_.find(it->first) != mass_.end()) nmass -= it->second;
    if (length_.find(it->first) != length_.end()) nlength -= it->second;
    if (volume_.find(it->first) != volume_.end()) nlength -= 3 * it->second;
    if (concentration_.find(it->first) != concentration_.end()) nconcentration -= it->second;
    if (amount_.find(it->first) != amount_.end()) namount -= it->second;
    if (temperature_.find(it->first) != temperature_.end()) ntemperature -= it->second;
  }

  return !(ntime || nmass || nlength || namount || nconcentration || ntemperature);
}


/* ******************************************************************
* Do two unit strings desribe the same physical quantity?
****************************************************************** */
bool
Units::CompareUnits(const std::string& unit1, const std::string& unit2)
{
  bool flag1, flag2;
  AtomicUnitForm auf1 = ComputeAtomicUnitForm_(unit1, &flag1);
  AtomicUnitForm auf2 = ComputeAtomicUnitForm_(unit2, &flag2);

  for (std::map<std::string, AtomicUnitForm>::iterator it = derived_.begin(); it != derived_.end();
       ++it) {
    auf1.replace(it->first, it->second);
    auf2.replace(it->first, it->second);
  }

  if (!flag1 || !flag2) return false;
  return CompareAtomicUnitForms_(auf1, auf2);
}


/* ******************************************************************
* Fancy output of time given in second
****************************************************************** */
std::string
Units::OutputTime(double val)
{
  double out, dmin, dtry, tmp;
  std::string unit("s");

  out = val;
  if (val == 0) return std::string("0 s");

  if (val > 0.0) {
    dmin = fabs(log10(val));

    for (auto it = time_.begin(); it != time_.end(); ++it) {
      tmp = val / it->second;
      dtry = fabs(log10(tmp));
      if (dtry < dmin) {
        dmin = dtry;
        out = tmp;
        unit = it->first;
      }
    }
  }

  std::stringstream ss;
  ss << out << " " << unit;
  return ss.str();
}


/* ******************************************************************
* Fancy output of mass given in kilograms
****************************************************************** */
std::string
Units::OutputMass(double val) const
{
  double out, dmin, dtry, tmp;
  std::string unit("kg");

  out = val;
  if (val == 0) return std::string("0 kg");

  if (val > 0.0) {
    dmin = fabs(log10(val));

    for (auto it = mass_.begin(); it != mass_.end(); ++it) {
      tmp = val / it->second;
      dtry = fabs(log10(tmp));
      if (dtry < dmin) {
        dmin = dtry;
        out = tmp;
        unit = it->first;
      }
    }
  }

  std::stringstream ss;
  ss << out << " " << unit;
  return ss.str();
}


/* ******************************************************************
* Fancy output of length given in meters
****************************************************************** */
std::string
Units::OutputLength(double val) const
{
  double out, dmin, dtry, tmp;
  std::string unit("m");

  out = val;
  if (val == 0) return std::string("0 m");

  if (val > 0.0) {
    dmin = fabs(log10(val));

    for (auto it = length_.begin(); it != length_.end(); ++it) {
      tmp = val / it->second;
      dtry = fabs(log10(tmp));
      if (dtry < dmin) {
        dmin = dtry;
        out = tmp;
        unit = it->first;
      }
    }
  }

  std::stringstream ss;
  ss << out << " " << unit;
  return ss.str();
}


/* ******************************************************************
* Fancy output of concentration
****************************************************************** */
std::string
Units::OutputConcentration(double val)
{
  std::string unit;
  unit = (system_.concentration == "molar") ? "mol/L" : "mol/m^3";

  std::stringstream ss;
  ss << val << " " << unit;
  return ss.str();
}

} // namespace Utils
} // namespace Amanzi
