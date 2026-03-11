/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Thermal conductivity for salt.
*/

#include "dbc.hh"
#include "errors.hh"

#include "NaCl_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor takes a parameter list to override defaulr values. This
* may be useful for unit tests.
******************************************************************* */
NaCl_ThermalConductivity::NaCl_ThermalConductivity(Teuchos::ParameterList& plist)
  : EOS_ThermalConductivity(plist), kref_(5.4), Tref_(300.00)
{
  kref_ = plist_.get<double>("reference conductivity", kref_);
  include_phi_ = plist_.get<bool>("include porosity dependence", true);
  clipping_ = true;
}


/* *******************************************************************
* Main functions
******************************************************************* */
double
NaCl_ThermalConductivity::ThermalConductivity(double p, double T, double phi)
{
  double Ts = Tref_ / T;
  double coef1 = kref_ * std::pow(Ts, 1.14);

  if (include_phi_) {
    double phi2 = phi * phi;
    double phi3 = phi2 * phi;
    double phi4 = phi3 * phi;
    double coef2 = 5.0 + 1.5 * phi - 136.0 * phi2 + 370.0 * phi3 - 270 * phi4;

    if (clipping_ && phi > 0.38815868) coef2 = 0.600967458596102;
    return coef1 * coef2 / 5.0;
  } else {
    return coef1;
  }
}


/* *******************************************************************
* Derivative wrt temperature
******************************************************************* */
double
NaCl_ThermalConductivity::DThermalConductivityDT(double p, double T, double phi)
{
  double Ts = Tref_ / T;
  double coef1 = -kref_ * 1.14 * std::pow(Ts, 0.14) / T / T;

  if (include_phi_) {
    double phi2 = phi * phi;
    double phi3 = phi2 * phi;
    double phi4 = phi3 * phi;
    double coef2 = 5.0 + 1.5 * phi - 136.0 * phi2 + 370.0 * phi3 - 270 * phi4;

    return coef1 * coef2 / 5.0;
  } else {
    return coef1;
  }
}


/* *******************************************************************
* Derivative wrt porosity
******************************************************************* */
double
NaCl_ThermalConductivity::DThermalConductivityDPhi(double p, double T, double phi)
{
  if (include_phi_) {
    double Ts = Tref_ / T;
    double coef1 = kref_ * std::pow(Ts, 1.14);

    double phi2 = phi * phi;
    double phi3 = phi2 * phi;
    double coef2 = 1.5 - 272.0 * phi + 1110.0 * phi2 - 1080 * phi3;

    if (clipping_ && phi > 0.38815868) coef2 = 0.0;

    return coef1 * coef2 / 5.0;
  } else {
    return 0.0;
  }
}

} // namespace AmanziEOS
} // namespace Amanzi
