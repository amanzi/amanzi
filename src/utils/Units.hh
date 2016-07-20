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

extern bool concentration_mol_liter;

// units for liter
typedef boost::units::make_scaled_unit<
    boost::units::si::volume, boost::units::scale<10, boost::units::static_rational<-3> > >::type liter;

// units for concentration
typedef boost::units::derived_dimension<
    boost::units::amount_base_dimension, 1,
    boost::units::length_base_dimension, -3>::type concentration_dimension;

typedef boost::units::unit<
    concentration_dimension, boost::units::si::system> concentration;

typedef boost::units::make_scaled_unit<
    concentration, boost::units::scale<10, boost::units::static_rational<3> > >::type concentration_amanzi;

// atomic representation of a derived unit
typedef std::map<std::string, int> AtomicUnitForm;

class Units {
 public:
  Units()
    : time_unit_("s"), 
      length_unit_("m"),
      mass_unit_("kg"),
      concentration_unit_("molar") {};
  Units(const std::string& concentration_unit) 
    : time_unit_("s"), 
      length_unit_("m"),
      mass_unit_("kg"),
      concentration_unit_(concentration_unit) { Init(); }
  ~Units() {};

  // conversion factors
  double pressure_factor() { return pressure_factor_; } 
  double concentration_factor() { return tcc_factor_; } 

  // output
  std::string print_tcc(double val, const std::string& system = "amanzi") {
    boost::units::quantity<concentration_amanzi> qval = val * concentration_amanzi();
    boost::units::quantity<concentration> qval_si(qval);
    std::stringstream ss;
    (system == "si") ? ss << qval_si : ss << qval;
    return ss.str();
  }

  void Init();

  void Init(Teuchos::ParameterList& plist) {
    concentration_unit_ = plist.get<std::string>("concentration", "molar");
    Init();
  }

  // conversion of units
  // -- deprecated
  double ConvertTime(double val, const std::string& unit); 

  // -- new
  double ConvertTime(double val, const std::string& in_unit,
                     const std::string& out_unit, bool& flag);

  double ConvertLength(double val, const std::string& in_unit,
                       const std::string& out_unit, bool& flag);

  double ConvertConcentration(double val, const std::string& in_unit,
                              const std::string& out_unit, double mol_mass, bool& flag);

  double ConvertDerivedUnit(double val, const std::string& in_unit,
                            const std::string& out_unit, double mol_mass, bool& flag);

  // access
  std::string time_unit() { return time_unit_; }
  std::string length_unit() { return length_unit_; }
  std::string mass_unit() { return mass_unit_; }
  std::string concentration_unit() { return concentration_unit_; }

 private:
  AtomicUnitForm ComputeAtomicUnitForm_(const std::string& unit, bool* flag);

 private:
  double pressure_factor_, tcc_factor_;

  std::map<std::string, boost::units::quantity<boost::units::si::time> > time_;
  std::map<std::string, boost::units::quantity<boost::units::si::length> > length_;
  std::map<std::string, boost::units::quantity<concentration> > concentration_;

  // default units
  std::string time_unit_;
  std::string length_unit_;
  std::string mass_unit_;
  std::string concentration_unit_;
};

}  // namespace Utils
}  // namespace Amanzi


namespace boost {
namespace units {

std::string symbol_string(const Amanzi::Utils::concentration_amanzi&);
std::string name_string(const Amanzi::Utils::concentration_amanzi&);

}  // namespace units
}  // namespace boost

#endif
