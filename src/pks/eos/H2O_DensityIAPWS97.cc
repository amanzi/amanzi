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
  Properties of Water and Steam (IAPWS), Industrial Formulation 1997.
*/

#include "H2O_DensityIAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_DensityIAPWS97::H2O_DensityIAPWS97(Teuchos::ParameterList& plist)
  : EOS_Density(plist) {
  eos_ = Teuchos::rcp(new IAPWS97(plist));
};


double
H2O_DensityIAPWS97::Density(double T, double p)
{
  double pMPa = p * 1e-6;
  return (eos_->ThermodynamicsPT(pMPa, T)).rho;
}


double
H2O_DensityIAPWS97::DDensityDT(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = eos_->ThermodynamicsPT(pMPa, T);
  if (prop.rgn == 3)
    return -prop.ap / (prop.v * prop.v * prop.bp);  // Helmholtz
  else
    return -prop.av / prop.v;  // Gibbs
}


double
H2O_DensityIAPWS97::DDensityDp(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = eos_->ThermodynamicsPT(pMPa, T);
  if (prop.rgn == 3)
    return 1.0e-6 / (prop.v * prop.v * pMPa * prop.bp);
  else
    return 1.0e-6 * prop.kt / prop.v;
}

} // namespace AmanziEOS
} // namespace Amanzi
