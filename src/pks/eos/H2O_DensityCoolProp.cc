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

#include "H2O_DensityCoolProp.hh"


namespace Amanzi {
namespace AmanziEOS {

using namespace CoolProp;

H2O_DensityCoolProp::H2O_DensityCoolProp(Teuchos::ParameterList& plist)
  : EOS_Density(plist)
{ 
  state_ = CoolProp::AbstractState::factory("HEOS", "Water");
}


double
H2O_DensityCoolProp::Density(double T, double p)
{
  state_->update(input_pairs::PT_INPUTS, p, T);
  return state_->rhomass();
};


double
H2O_DensityCoolProp::DDensityDT(double T, double p)
{
  state_->update(input_pairs::PT_INPUTS, p, T);
  return state_->first_partial_deriv(iDmass, iT, iP);
};


double
H2O_DensityCoolProp::DDensityDp(double T, double p)
{
  state_->update(input_pairs::PT_INPUTS, p, T);
  return state_->first_partial_deriv(iDmass, iP, iT);
};


CoolProp::phases
H2O_DensityCoolProp::get_phase(double T, double p) {
  state_->update(input_pairs::PT_INPUTS, p, T);
  return state_->phase();
}

} // namespace AmanziEOS
} // namespace Amanzi
