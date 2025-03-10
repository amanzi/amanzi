/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Fugacity for liquid water is defined as the saturated water vapor pressure.

*/

#ifndef AMANZI_FUGACITY_SATURATED_VAPOR_HH_
#define AMANZI_FUGACITY_SATURATED_VAPOR_HH_

#include <cmath>

#include "dbc.hh"

#include "Fugacity.hh"

namespace Amanzi {
namespace Multiphase {

class Fugacity_SaturatedVapor : public Fugacity {
 public:
  Fugacity_SaturatedVapor(const Teuchos::ParameterList& plist) 
   : ka0_(16.635764),
     ka_(-6096.9385),
     kb_(-2.7111933e-2),
     kc_(1.673952e-5),
     kd_(2.433502){};

  ~Fugacity_SaturatedVapor() {};

  virtual double Value(double T) override {
    if (T < 100.0 || T > 373.15) AMANZI_ASSERT(false);
    return 100.0 * std::exp(ka0_ + ka_ / T + (kb_ + kc_ * T) * T + kd_ * std::log(T));
  }

 private:
   double ka0_, ka_, kb_, kc_, kd_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
