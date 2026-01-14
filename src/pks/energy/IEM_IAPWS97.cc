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

#include "CommonDefs.hh"
#include "IEM_IAPWS97.hh"

namespace Amanzi {
namespace Energy {

static constexpr double cfactor = 1000.0 * CommonDefs::MOLAR_MASS_H2O;

IEM_IAPWS97::IEM_IAPWS97(Teuchos::ParameterList& plist)
{
  eos_ = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));
}


double
IEM_IAPWS97::InternalEnergy(double T, double p)
{
  double pMPa = p * 1e-6;
  return cfactor * eos_->ThermodynamicsPT(pMPa, T).u;
}


double
IEM_IAPWS97::DInternalEnergyDT(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = eos_->ThermodynamicsPT(pMPa, T);
  if (prop.rgn == 3)
    return cfactor * (pMPa * (T * prop.ap - 1.0) - prop.cv) / pMPa / prop.bp;
  else
    return -cfactor * prop.v * (pMPa * prop.kt - T * prop.av);
}


double
IEM_IAPWS97::DInternalEnergyDp(double T, double p)
{
  double pMPa = p * 1e-6;
  auto prop = eos_->ThermodynamicsPT(pMPa, T);
  if (prop.rgn == 3)
    return cfactor * (1.0 - T * prop.ap) / prop.bp;
  else
    return cfactor * prop.cp - pMPa * prop.v * prop.av;
}

} // namespace Energy
} // namespace Amanzi
