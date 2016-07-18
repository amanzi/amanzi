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

class Units {
 public:
  Units() {};
  Units(const std::string& conc_units) { Init(conc_units); }
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

  void Init(const std::string& conc_units);

  void Init(Teuchos::ParameterList& plist) {
    Init(plist.get<std::string>("concentration", "molar"));
  }

 // conversion of units
 double convert_time(double t, const std::string& unit); 

 public:
  std::map<std::string, boost::units::quantity<boost::units::si::length> > length_;

 private:
  double pressure_factor_, tcc_factor_;
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
