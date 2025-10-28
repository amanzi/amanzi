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
  state_->update(input_pairs::PT_INPUTS, p, T);
  return state_->viscosity();
};


double
H2O_ViscosityCoolProp::DViscosityDT(double T, double p)
{
  state_->update(input_pairs::PT_INPUTS, p, T);
  double mu = state_->viscosity();

  double dT = 0.01;
  state_->update(PT_INPUTS, p, T + dT);
  return (state_->viscosity() - mu) / dT;

};


double
H2O_ViscosityCoolProp::DViscosityDp(double T, double p)
{
  state_->update(input_pairs::PT_INPUTS, p, T);
  double mu = state_->viscosity();

  double dp = 1000.0;
  state_->update(PT_INPUTS, p + dp, T);
  return (state_->viscosity() - mu) / dp;
};


CoolProp::phases
H2O_ViscosityCoolProp::get_phase(double T, double p) {
  state_->update(input_pairs::PT_INPUTS, p, T);
  return state_->phase();
}

} // namespace AmanziEOS
} // namespace Amanzi
