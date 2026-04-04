/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  EOS for liquid water based on the International Association for the 
  Properties of Water and Steam (IAPWS) formulation 1995.
*/

#include "H2O_DensityIAPWS95.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_DensityIAPWS95::H2O_DensityIAPWS95(Teuchos::ParameterList& plist)
  : EOS_Density(plist) {
  eos_ = Teuchos::rcp(new IAPWS95(plist));
};


double
H2O_DensityIAPWS95::Density(double T, double p)
{
  double pMPa = p * 1e-6;
  return std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).rho;
}


double
H2O_DensityIAPWS95::DDensityDT(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = std::get<0>(eos_->ThermodynamicsPT(pMPa, T));
  return -prop.ap / (prop.v * prop.v * prop.bp);  // Helmholtz
}


double
H2O_DensityIAPWS95::DDensityDp(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = std::get<0>(eos_->ThermodynamicsPT(pMPa, T));
  return 1.0e-6 / (prop.v * prop.v * pMPa * prop.bp);
}

} // namespace AmanziEOS
} // namespace Amanzi
