/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Viscosity for liquid water consistent with IAPWS97 formulation.
  http://www.iapws.org/relguide/viscosity.html, formulas (10)-(12).
*/

#include "errors.hh"
#include "H2O_ViscosityIAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_ViscosityIAPWS97::H2O_ViscosityIAPWS97(Teuchos::ParameterList& plist)
  : EOS_Viscosity(plist)
{
  eos_ = Teuchos::rcp(new IAPWS97(plist));
}


double
H2O_ViscosityIAPWS97::Viscosity(double T, double p)
{
  double pMPa = p * 1e-6;
  return eos_->ThermodynamicsPT(pMPa, T).mu;
}


double
H2O_ViscosityIAPWS97::DViscosityDT(double T, double p)
{
  double pMPa = p * 1e-6;
  double mu1 = eos_->ThermodynamicsPT(pMPa, T).mu;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dT = eps * T;
  double mu2 = eos_->ThermodynamicsPT(pMPa, T + dT).mu;

  return (mu2 - mu1) / dT;
};


double
H2O_ViscosityIAPWS97::DViscosityDp(double T, double p)
{
  double pMPa = p * 1e-6;
  double mu1 = eos_->ThermodynamicsPT(pMPa, T).mu;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dp = eps * pMPa;
  double mu2 = eos_->ThermodynamicsPT(pMPa + dp, T).mu;

  return (mu2 - mu1) / dp * 1e-6;
};

} // namespace AmanziEOS
} // namespace Amanzi
