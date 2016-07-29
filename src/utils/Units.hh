/*
  Utils

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#ifndef AMANZI_UNITS_HH_
#define AMANZI_UNITS_HH_

#include <iostream>
#include <cstdio>
#include <iomanip>

#include <boost/units/scaled_base_unit.hpp>
#include <boost/units/derived_dimension.hpp>
#include <boost/units/io.hpp>
#include <boost/units/make_scaled_unit.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/unit.hpp>

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Utils {

typedef boost::units::derived_dimension<
    boost::units::amount_base_dimension, 1,
    boost::units::length_base_dimension, -3>::type concentration_dimension;

typedef boost::units::unit<
    concentration_dimension, boost::units::si::system> concentration;


// ------------------------------------------------------------------
// Auxiliary class: atomic representation of a derived unit. 
// It may contain other derived units.
// ------------------------------------------------------------------
typedef std::map<std::string, int> UnitData;

class AtomicUnitForm {
 public:
  AtomicUnitForm() {};
  AtomicUnitForm(const std::string& key1, int i1) {
    data_[key1] = i1;
  }
  AtomicUnitForm(const std::string& key1, int i1,
                 const std::string& key2, int i2) {
    data_[key1] = i1;
    data_[key2] = i2;
  }
  AtomicUnitForm(const std::string& key1, int i1,
                 const std::string& key2, int i2,
                 const std::string& key3, int i3) {
    data_[key1] = i1;
    data_[key2] = i2;
    data_[key3] = i3;
  }
  ~AtomicUnitForm() {};

  // elementary operations with reduced representations
  // -- optional replacement of a key
  AtomicUnitForm& replace(const std::string& key, const AtomicUnitForm& aut) {
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
  void print() {
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
              const std::string& concentration_) :
      time(time_),
      mass(mass_),
      length(length_),
      concentration(concentration_) {};

  std::string time;
  std::string mass;
  std::string length;
  std::string concentration;
};


class Units {
 public:
  Units() : system_("s", "kg", "m", "molar") { Init(); }
  Units(const std::string& concentration_unit) : system_("s", "kg", "m", "molar") { 
    system_.concentration = concentration_unit;
    Init();
  }
  ~Units() {};

  // main members
  void Init();

  void Init(Teuchos::ParameterList& plist) {
    system_.concentration = plist.get<std::string>("concentration", "molar");
    Init();
  }

  // conversion factors
  double concentration_factor() { return concentration_factor_; } 

  // conversion of units
  // -- data
  double ConvertTime(double val, const std::string& in_unit,
                     const std::string& out_unit, bool& flag);

  double ConvertMass(double val, const std::string& in_unit,
                     const std::string& out_unit, bool& flag);

  double ConvertLength(double val, const std::string& in_unit,
                       const std::string& out_unit, bool& flag);

  double ConvertConcentration(double val, const std::string& in_unit,
                              const std::string& out_unit, double mol_mass, bool& flag);

  double ConvertUnitD(double val, const std::string& in_unit,
                      const std::string& out_unit, double mol_mass, bool& flag);

  // -- strings
  std::string ConvertUnitS(const std::string& in_unit, const UnitsSystem& system);

  // -- fancy output
  std::string OutputTime(double val);
  std::string OutputConcentration(double val);

  // access
  UnitsSystem& system() { return system_; }

 private:
  AtomicUnitForm ComputeAtomicUnitForm_(const std::string& unit, bool* flag);

 private:
  double concentration_factor_;

  std::map<std::string, boost::units::quantity<boost::units::si::time> > time_;
  std::map<std::string, boost::units::quantity<boost::units::si::mass> > mass_;
  std::map<std::string, boost::units::quantity<boost::units::si::length> > length_;
  std::map<std::string, boost::units::quantity<boost::units::si::volume> > volume_;
  std::map<std::string, boost::units::quantity<concentration> > concentration_;

  std::map<std::string, AtomicUnitForm> derived_;

  // default Amanzi's units
  UnitsSystem system_;
};

}  // namespace Utils
}  // namespace Amanzi

#endif
