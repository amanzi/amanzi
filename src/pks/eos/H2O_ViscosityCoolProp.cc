/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  EOS for liquid water from package CoolProp
*/

#include <CoolProp.h>

#include "H2O_ViscosityCoolProp.hh"


namespace Amanzi {
namespace AmanziEOS {

using namespace CoolProp;

H2O_ViscosityCoolProp::H2O_ViscosityCoolProp(Teuchos::ParameterList& plist)
  : EOS_Viscosity(plist)
{ 
  state_ = CoolProp::AbstractState::factory("HEOS", "Water");
}


double
H2O_ViscosityCoolProp::Viscosity(double T, double p)
{
  if (update(T, p)) return state_->viscosity();
  return 0.0;
};


double
H2O_ViscosityCoolProp::DViscosityDT(double T, double p)
{
  double mu(0.0);
  if (update(T, p)) mu = state_->viscosity();

  double dT = 0.01;
  if (update(T + dT, p)) return (state_->viscosity() - mu) / dT;
  return 0.0;
};


double
H2O_ViscosityCoolProp::DViscosityDp(double T, double p)
{
  double mu(0.0);
  if (update(T, p)) mu = state_->viscosity();

  double dp = 1000.0;
  if (update(T, p + dp)) return (state_->viscosity() - mu) / dp;
  return 0.0;
};


CoolProp::phases
H2O_ViscosityCoolProp::get_phase(double T, double p) {
  if (update(T, p)) return state_->phase();
  return CoolProp::phases::iphase_liquid;
}


bool
H2O_ViscosityCoolProp::update(double T, double p)
{
  try {
    state_->update(input_pairs::PT_INPUTS, p, T);  
    return true;
  } 
  catch (const CoolProp::ValueError &e) {
    ierr_ = 1;
    std::stringstream ss;
    ss << e.what();
    error_msg_ = ss.str();
  }
  catch (const std::exception &e) {
    ierr_ = 1;
    std::stringstream ss;
    ss << e.what();
    error_msg_ = ss.str();
  }
  return false;
}

} // namespace AmanziEOS
} // namespace Amanzi
