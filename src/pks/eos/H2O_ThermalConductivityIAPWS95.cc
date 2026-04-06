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
  Properties of Water and Steam (IAPWS) Formulation 1995.
  http://www.iapws.org/relguide/ThCond.html
*/

#include "dbc.hh"
#include "errors.hh"

#include "H2O_ThermalConductivityIAPWS95.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor sets up EOS
******************************************************************* */
H2O_ThermalConductivityIAPWS95::H2O_ThermalConductivityIAPWS95(Teuchos::ParameterList& plist)
  : EOS_ThermalConductivity(plist)
{
  eos_ = Teuchos::rcp(new IAPWS95(plist));
}


/* *******************************************************************
* Safe-guarded state evaluation for parallel execution.  
******************************************************************* */
double
H2O_ThermalConductivityIAPWS95::ThermalConductivity(double p, double T, double phi)
{
  double k, pMPa(p * 1e-6);
  try {
    k = std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).k;
  } catch (...) {
    ierr_ = 1;
    std::stringstream ss;
    ss << "invalid T=" << T << " conductivity=" << k;
    error_msg_ = ss.str();
  }
  return k;
}


/* *******************************************************************
* Derivative is based on finite difefreces
******************************************************************* */
double
H2O_ThermalConductivityIAPWS95::DThermalConductivityDp(double p, double T, double phi)
{
  double pMPa = p * 1e-6;
  double k1 = std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).k;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dp = eps * pMPa;
  double k2 = std::get<0>(eos_->ThermodynamicsPT(pMPa + dp, T)).k;

  return (k2 - k1) / dp * 1e-6;
}


/* *******************************************************************
* Derivative is based on finite difefreces
******************************************************************* */
double
H2O_ThermalConductivityIAPWS95::DThermalConductivityDT(double p, double T, double phi)
{
  double pMPa = p * 1e-6;
  double k1 = std::get<0>(eos_->ThermodynamicsPT(pMPa, T)).k;

  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double dT = eps * T;
  double k2 = std::get<0>(eos_->ThermodynamicsPT(pMPa, T + dT)).k;

  return (k2 - k1) / dT;
}

} // namespace AmanziEOS
} // namespace Amanzi
