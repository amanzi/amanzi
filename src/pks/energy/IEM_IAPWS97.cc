/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  EOS for liquid water based on the International Association for the 
  Properties of Water and Steam (IAPWS), Industrial Formulation 1997.
*/

#include "IEM_IAPWS97.hh"

namespace Amanzi {
namespace Energy {

IEM_IAPWS97::IEM_IAPWS97(Teuchos::ParameterList& plist)
{
  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));
}


double
IEM_IAPWS97::InternalEnergy(double T, double p)
{
  double pMPa = p * 1e-6;
  return 1000.0 * eos_->ThermodynamicsPT(pMPa, T).u;
}


double
IEM_IAPWS97::DInternalEnergyDT(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = eos_->ThermodynamicsPT(pMPa, T);
  if (prop.rgn == 3)
    return 1000.0 * (pMPa * (T * prop.ap - 1.0) - prop.cv) / pMPa / prop.bp;
  else
    return -1000.0 * prop.v * (pMPa * prop.kt - T * prop.av);
}


double
IEM_IAPWS97::DInternalEnergyDp(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = eos_->ThermodynamicsPT(pMPa, T);
  if (prop.rgn == 3)
    return 1000.0 * (1.0 - T * prop.ap) / prop.bp;
  else
    return 1000.0 * prop.cp - pMPa * prop.v * prop.av;
}

} // namespace Energy
} // namespace Amanzi
