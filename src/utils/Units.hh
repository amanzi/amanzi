/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*!

Amanzi's internal default units are SI units except for the concentration.

Curently format of units is very rigid, we need to use multiplication and
division operators, e.g. Pa*s and kg/s. However, thos allows to mix and match,
e.g. m/s^2 or ms^-2 or s^-2*m.

* `"concentration`" [string] defines units for concentration. Available options
  are `"molar`" (default) which is `"mol/L`" and `"SI`" which is `"mol/m^3`".

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="units">
    <Parameter name="length" type="string" value="m"/>
    <Parameter name="time" type="string" value="s"/>
    <Parameter name="mass" type="string" value="kg"/>
    <Parameter name="temperature" type="string" value="K"/>
    <Parameter name="concentration" type="string" value="molar"/>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_UNITS_HH_
#define AMANZI_UNITS_HH_

#include <iostream>
#include <cstdio>
#include <iomanip>

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Utils {

// ------------------------------------------------------------------
// Auxiliary class: atomic representation of a derived unit.
// It may contain other derived units.
// ------------------------------------------------------------------
typedef std::map<std::string, int> UnitData;

class AtomicUnitForm {
 public:
  AtomicUnitForm(){};
  AtomicUnitForm(const std::string& key1, int i1) { data_[key1] = i1; }
  AtomicUnitForm(const std::string& key1, int i1, const std::string& key2, int i2)
  {
    data_[key1] = i1;
    data_[key2] = i2;
  }
  AtomicUnitForm(const std::string& key1,
                 int i1,
                 const std::string& key2,
                 int i2,
                 const std::string& key3,
                 int i3)
  {
    data_[key1] = i1;
    data_[key2] = i2;
    data_[key3] = i3;
  }
  ~AtomicUnitForm(){};

  // elementary operations with reduced representations
  // -- optional replacement of a key
  AtomicUnitForm& replace(const std::string& key, const AtomicUnitForm& aut)
  {
    const UnitData& data = aut.data();

    UnitData::iterator it1;
    UnitData::const_iterator it2;

    if ((it1 = data_.find(key)) != data_.end()) {
      int i = it1->second;
      data_.erase(it1);

      for (it2 = data.begin(); it2 != data.end(); ++it2) {
        if ((it1 = data_.find(it2->first)) != data_.end()) {
          it1->second += i * it2->second;
        } else {
          data_[it2->first] = i * it2->second;
        }
      }
    }
    return *this;
  }

  // access
  const UnitData& data() const { return data_; }
  UnitData& data() { return data_; }

  // development
  void print()
  {
    std::cout << "AtomicUnitForm:" << std::endl;
    for (UnitData::iterator it = data_.begin(); it != data_.end(); ++it) {
      std::cout << "  " << it->first << " " << it->second << std::endl;
    }
  }

 private:
  UnitData data_;
};


// ------------------------------------------------------------------
// Base class for converting units.
// ------------------------------------------------------------------
struct UnitsSystem {
  UnitsSystem(const std::string& time_,
              const std::string& mass_,
              const std::string& length_,
              const std::string& concentration_,
              const std::string& amount_,
              const std::string& temperature_)
    : time(time_),
      mass(mass_),
      length(length_),
      concentration(concentration_),
      amount(amount_),
      temperature(temperature_){};

  std::string time;
  std::string mass;
  std::string length;
  std::string concentration;
  std::string amount;
  std::string temperature;
};


class Units {
 public:
  Units() : system_("s", "kg", "m", "molar", "mol", "K") { Init(); }
  Units(const std::string& concentration_unit) : system_("s", "kg", "m", "molar", "mol", "K")
  {
    system_.concentration = concentration_unit;
    Init();
  }
  ~Units(){};

  // main members
  void Init();

  void Init(Teuchos::ParameterList& plist)
  {
    system_.time = plist.get<std::string>("time", "s");
    system_.mass = plist.get<std::string>("mass", "kg");
    system_.length = plist.get<std::string>("length", "m");
    system_.concentration = plist.get<std::string>("concentration", "molar");
    system_.temperature = plist.get<std::string>("temperature", "K");
    Init();
  }

  // conversion factors
  double concentration_factor() { return concentration_factor_; }

  // conversion/comparison of units
  // -- data
  double
  ConvertTime(double val, const std::string& in_unit, const std::string& out_unit, bool& flag);

  double
  ConvertMass(double val, const std::string& in_unit, const std::string& out_unit, bool& flag);

  double
  ConvertLength(double val, const std::string& in_unit, const std::string& out_unit, bool& flag);

  double ConvertConcentration(double val,
                              const std::string& in_unit,
                              const std::string& out_unit,
                              double mol_mass,
                              bool& flag);

  double ConvertTemperature(double val,
                            const std::string& in_unit,
                            const std::string& out_unit,
                            bool& flag);

  double ConvertUnitD(double val,
                      const std::string& in_unit,
                      const std::string& out_unit,
                      double mol_mass,
                      bool& flag);

  // -- strings
  std::string ConvertUnitS(const std::string& in_unit, const UnitsSystem& system);
  bool CompareUnits(const std::string& unit1, const std::string& unit2);
  std::string MultiplyUnits(const std::string& unit1, const std::string& unit2);
  std::string DivideUnits(const std::string& unit1, const std::string& unit2);

  // fancy output
  std::string OutputTime(double val);
  std::string OutputMass(double val) const;
  std::string OutputLength(double val) const;
  std::string OutputConcentration(double val);

  // error checking
  bool IsValidTime(const std::string& unit) const { return time_.find(unit) != time_.end(); }
  std::string ValidTimeStrings() const
  {
    std::stringstream valids;
    for (auto v : time_) valids << "\"" << v.first << "\",";
    auto valids_str = valids.str();
    valids_str.pop_back(); // remove the last comma
    return valids_str;
  }

  bool IsValidMass(const std::string& unit) const { return mass_.find(unit) != mass_.end(); }
  std::string ValidMassStrings() const
  {
    std::stringstream valids;
    for (auto v : mass_) valids << "\"" << v.first << "\",";
    auto valids_str = valids.str();
    valids_str.pop_back(); // remove the last comma
    return valids_str;
  }

  bool IsValidLength(const std::string& unit) const { return length_.find(unit) != length_.end(); }
  std::string ValidLengthStrings() const
  {
    std::stringstream valids;
    for (auto v : length_) valids << "\"" << v.first << "\",";
    auto valids_str = valids.str();
    valids_str.pop_back(); // remove the last comma
    return valids_str;
  }

  bool IsValidConcentration(const std::string& unit) const
  {
    return concentration_.find(unit) != concentration_.end();
  }
  std::string ValidConcentrationStrings() const
  {
    std::stringstream valids;
    for (auto v : concentration_) valids << "\"" << v.first << "\",";
    auto valids_str = valids.str();
    valids_str.pop_back(); // remove the last comma
    return valids_str;
  }

  bool IsValidTemperature(const std::string& unit) const
  {
    return temperature_.find(unit) != temperature_.end();
  }
  std::string ValidTemperatureStrings() const
  {
    std::stringstream valids;
    for (auto v : temperature_) valids << "\"" << v.first << "\",";
    auto valids_str = valids.str();
    valids_str.pop_back(); // remove the last comma
    return valids_str;
  }

  // access
  UnitsSystem& system() { return system_; }

 private:
  AtomicUnitForm StringToAtomicUnitForm_(const std::string& unit, bool* flag);
  std::string AtomicUnitFormToString_(const AtomicUnitForm& auf);

  bool CompareAtomicUnitForms_(const AtomicUnitForm& auf1, const AtomicUnitForm& auf2);
  AtomicUnitForm MultiplyAtomicUnitForms_(const AtomicUnitForm& auf1, const AtomicUnitForm& auf2);
  AtomicUnitForm DivideAtomicUnitForms_(const AtomicUnitForm& auf1, const AtomicUnitForm& auf2);

 private:
  double concentration_factor_;

  std::map<std::string, double> time_;
  std::map<std::string, double> mass_;
  std::map<std::string, double> length_;
  std::map<std::string, double> volume_;
  std::map<std::string, double> concentration_;
  std::map<std::string, double> amount_;
  std::map<std::string, double> temperature_;

  std::map<std::string, AtomicUnitForm> derived_;

  // default Amanzi's units
  UnitsSystem system_;
};


} // namespace Utils
} // namespace Amanzi

#endif
