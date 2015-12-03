/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#ifndef AMANZI_UNITS_HH_
#define AMANZI_UNITS_HH_

#include <boost/units/make_scaled_unit.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/unit.hpp>
#include <boost/units/derived_dimension.hpp>

#include "Teuchos_ParameterList.hpp"

namespace Utils {

// units for Amanzi concentration function
typedef boost::units::derived_dimension<
    boost::units::amount_base_dimension, 1,
    boost::units::length_base_dimension, -3>::type concentration_si;

typedef boost::units::make_scaled_unit<
    boost::units::si::volume, boost::units::scale<10, boost::units::static_rational<-3> > >::type liter;

class Units {
 public:
  Units() {};
  Units(const std::string& conc_units) { Init(conc_units); }
  ~Units() {};

  // conversion factors
  inline double pressure_factor() { return pressure_factor_; } 
  inline double concentration_factor() { return tcc_factor_; } 

  // output
  /*
  friend std::ostream& operator << (std::ostream& os, double val) {
    os << (val * tcc_si_);
    return os;
  }
  */

  void Init(const std::string& conc_units) {
    pressure_factor_ = conversion_factor(boost::units::si::pressure(), pressure_si_);
    if (conc_units == "molar") {
      tcc_factor_ = conversion_factor(boost::units::si::amount() / liter(), tcc_si_);
    } else {
      tcc_factor_ = 1.0;
    }
  }

  inline void Init(Teuchos::ParameterList& plist) {
    Init(plist.get<std::string>("concentration", "molar"));
  }

 private:
  boost::units::si::pressure pressure_si_;
  boost::units::unit<concentration_si, boost::units::si::system> tcc_si_;

  double pressure_factor_, tcc_factor_;
};

}  // namespace Utils

#endif
