/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Thermal conductivity for ordinary liquid water.
*/

#include "dbc.hh"
#include "errors.hh"

#include "H2O_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor takes a parameter list to override defaulr values. This
* may be useful for unit tests.
******************************************************************* */
H2O_ThermalConductivity::H2O_ThermalConductivity(Teuchos::ParameterList& plist)
  : EOS_ThermalConductivity(plist),
    ka0_(-1.48445),
    ka1_(4.12292),
    ka2_(-1.63866),
    kref_(0.6065),
    Tref_(298.15)
{
  kref_ = plist_.get<double>("reference conductivity", kref_);
  Tref_ = plist_.get<double>("reference temperature", Tref_);

  if (plist_.isParameter("polynomial expansion")) {
    Teuchos::Array<double> kai = plist_.get<Teuchos::Array<double>>("polynomial expansion");
    AMANZI_ASSERT(kai.size() == 3);

    ka0_ = kai[0];
    ka1_ = kai[1];
    ka2_ = kai[2];
  }

  ka0_ *= kref_;
  ka1_ *= kref_;
  ka2_ *= kref_;
}


/* *******************************************************************
* Main routine uses pre-computed coefficients.
******************************************************************* */
double
H2O_ThermalConductivity::ThermalConductivity(double T, double phi)
{
  double Ts = T / Tref_;
  double k = ka0_ + (ka1_ + ka2_ * Ts) * Ts;

  ierr_ = 0;
  if (k < 1.0e-16) {
    ierr_ = 1;
    std::stringstream ss;
    ss << "invalid T=" << T << " conductivity=" << k;
    error_msg_ = ss.str();
  }
  return k;
}


/* *******************************************************************
* Main routine uses pre-computed coefficients.
******************************************************************* */
double
H2O_ThermalConductivity::DThermalConductivityDT(double T, double phi)
{
  double Ts = T / Tref_;
  double dkdT = (ka1_ + 2 * ka2_ * Ts) / Tref_;

  return dkdT;
}

} // namespace AmanziEOS
} // namespace Amanzi
