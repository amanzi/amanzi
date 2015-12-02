/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#ifndef AMANZI_UNITS_HH_
#define AMAZNI_UNITS_HH_

#include <boost/units/make_scaled_unit.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/unit.hpp>
#include <boost/units/derived_dimension.hpp>

#include "Teuchos_ParameterList.hpp"

namespace Utils {

// derived unit for concentration density functions
typedef boost::units::make_scaled_unit<
    boost::units::si::volume, boost::units::scale<10, boost::units::static_rational<-3> > >::type liter;

typedef boost::units::derived_dimension<boost::units::amount_base_dimension, 1,
                                        boost::units::length_base_dimension, -3>::type concentration_si;

class Units {
 public:
  Units() {};
  Units(Teuchos::ParameterList& plist) {};
  ~Units() {};

  double concentration_factor() {
    return 1.0;
    return conversion_factor(boost::units::si::amount() / liter(),
                             boost::units::si::amount() / boost::units::si::volume());
  } 

 private:
};

}  // namespace Utils

#endif
