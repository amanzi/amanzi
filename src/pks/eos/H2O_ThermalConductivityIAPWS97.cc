/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  EOS for liquid water based on the International Association for the 
  Properties of Water and Steam (IAPWS), Industrial Formulation 1997.
  http://www.iapws.org/relguide/ThCond.html
*/

#include "dbc.hh"
#include "errors.hh"

#include "H2O_ThermalConductivityIAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor sets up EOS
******************************************************************* */
H2O_ThermalConductivityIAPWS97::H2O_ThermalConductivityIAPWS97(Teuchos::ParameterList& plist)
  : EOS_ThermalConductivity(plist)
{
  eos_ = Teuchos::rcp(new IAPWS97(plist));
}


/* *******************************************************************
* 
******************************************************************* */
double
H2O_ThermalConductivityIAPWS97::ThermalConductivity(double p, double T, double phi)
{
  double pMPa = p * 1e-6;
  return eos_->ThermodynamicsPT(pMPa, T).k;
}


/* *******************************************************************
* Derivative is based on finite difefreces
******************************************************************* */
double
H2O_ThermalConductivityIAPWS97::DThermalConductivityDp(double p, double T, double phi)
{
  double pMPa = p * 1e-6;
  double k1 = eos_->ThermodynamicsPT(pMPa, T).k;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dp = eps * pMPa;
  double k2 = eos_->ThermodynamicsPT(pMPa + dp, T).k;

  return (k2 - k1) / dp * 1e-6;
}


/* *******************************************************************
* Derivative is based on finite difefreces
******************************************************************* */
double
H2O_ThermalConductivityIAPWS97::DThermalConductivityDT(double p, double T, double phi)
{
  double pMPa = p * 1e-6;
  double k1 = eos_->ThermodynamicsPT(pMPa, T).k;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dT = eps * T;
  double k2 = eos_->ThermodynamicsPT(pMPa, T + dT).k;

  return (k2 - k1) / dT;
}

} // namespace AmanziEOS
} // namespace Amanzi
