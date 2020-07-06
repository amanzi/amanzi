/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Thermal conductivity for liquid water.
*/

#include "dbc.hh"
#include "errors.hh"

#include "ThermalConductivity_Water.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor takes a parameter list to override defaulr values. This
* may be useful for unit tests.
******************************************************************* */
ThermalConductivity_Water::ThermalConductivity_Water(Teuchos::ParameterList& eos_plist) :
    eos_plist_(eos_plist),
    ka0_(-1.48445),
    ka1_(4.12292),
    ka2_(-1.63866),
    kref_(0.6065),
    Tref_(298.15)
{
  InitializeFromPlist_();

  ka0_ *= kref_;
  ka1_ *= kref_;
  ka2_ *= kref_;
}


/* *******************************************************************
* Main routine uses pre-computed coefficients.
******************************************************************* */
double ThermalConductivity_Water::ThermalConductivity(double T)
{
  double Ts = T / Tref_;
  double k = ka0_ + (ka1_ + ka2_ * Ts) * Ts;

  if (k < 1.e-16) {
    std::cout << "Invalid temperature, T=" << T << " conductivity=" << k << std::endl;
    Exceptions::amanzi_throw(Errors::CutTimeStep());
  }
  return k;
}


/* *******************************************************************
* Main routine uses pre-computed coefficients.
******************************************************************* */
double ThermalConductivity_Water::DThermalConductivityDT(double T)
{
  double Ts = T / Tref_;
  double dkdT = ka1_ + 2 * ka2_ * Ts;

  return dkdT;
}


/* *******************************************************************
* Here we can override some parameters.
******************************************************************* */
void ThermalConductivity_Water::InitializeFromPlist_()
{
  kref_ = eos_plist_.get<double>("thermal conductivity of liquid", kref_);
  Tref_ = eos_plist_.get<double>("reference temperature", Tref_);

  if (eos_plist_.isParameter("polynomial expansion")) {
    Teuchos::Array<double> kai = eos_plist_.get<Teuchos::Array<double> >("polynomial expansion");
    AMANZI_ASSERT(kai.size() == 3);

    ka0_ = kai[0]; 
    ka1_ = kai[1]; 
    ka2_ = kai[2]; 
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi
