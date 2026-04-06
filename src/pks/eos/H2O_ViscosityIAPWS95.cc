/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Viscosity for liquid water consistent with IAPWS95 formulation.
  http://www.iapws.org/relguide/viscosity.html, formulas (10)-(12).
*/

#include "errors.hh"
#include "H2O_ViscosityIAPWS95.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_ViscosityIAPWS95::H2O_ViscosityIAPWS95(Teuchos::ParameterList& plist)
  : EOS_Viscosity(plist)
{
  eos_ = Teuchos::rcp(new IAPWS95(plist));
}


double
H2O_ViscosityIAPWS95::Viscosity(double T, double p)
{
  double pMPa = p * 1e-6;
  return std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).mu;
}


double
H2O_ViscosityIAPWS95::DViscosityDT(double T, double p)
{
  double pMPa = p * 1e-6;
  double mu1 = std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).mu;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dT = eps * T;
  double mu2 = std::get<0>(eos_->ThermodynamicsPT(pMPa, T + dT)).mu;

  return (mu2 - mu1) / dT;
};


double
H2O_ViscosityIAPWS95::DViscosityDp(double T, double p)
{
  double pMPa = p * 1e-6;
  double mu1 = std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).mu;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dp = eps * pMPa;
  double mu2 = std::get<0>(eos_->ThermodynamicsPT(pMPa + dp, T)).mu;

  return (mu2 - mu1) / dp * 1e-6;
};

} // namespace AmanziEOS
} // namespace Amanzi
