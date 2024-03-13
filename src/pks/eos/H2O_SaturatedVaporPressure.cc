/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Saturated vapor pressure for vapor over water or ice, Sonntag (1990)
*/

#include <cmath>

#include "errors.hh"

#include "H2O_SaturatedVaporPressure.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_SaturatedVaporPressure::H2O_SaturatedVaporPressure(Teuchos::ParameterList& plist)
  : plist_(plist),
    ka0_(16.635764),
    ka_(-6096.9385),
    kb_(-2.7111933e-2),
    kc_(1.673952e-5),
    kd_(2.433502){};


double
H2O_SaturatedVaporPressure::Pressure(double T)
{
  ierr_ = 0;
  if (T < 100.0 || T > 373.0) {
    ierr_ = 1;
    std::stringstream ss;
    ss << "invalid T = " << T;
    error_msg_ = ss.str();
    return 0.0;
  } else {
    return 100.0 * exp(ka0_ + ka_ / T + (kb_ + kc_ * T) * T + kd_ * log(T));
  }
}


double
H2O_SaturatedVaporPressure::DPressureDT(double T)
{
  ierr_ = 0;
  if (T < 100.0 || T > 373.0) {
    ierr_ = 1;
    std::stringstream ss;
    ss << "invalid T = " << T;
    error_msg_ = ss.str();
    return 0.0;
  } else {
    return Pressure(T) * (-ka_ / (T * T) + kb_ + 2 * kc_ * T + kd_ / T);
  }
}

} // namespace AmanziEOS
} // namespace Amanzi
