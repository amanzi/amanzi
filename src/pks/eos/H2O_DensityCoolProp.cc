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
  if (update(T, p)) return state_->rhomass();
  return 0.0;
};


double
H2O_DensityCoolProp::DDensityDT(double T, double p)
{
  if (update(T, p)) return state_->first_partial_deriv(iDmass, iT, iP);
  return 0.0;
};


double
H2O_DensityCoolProp::DDensityDp(double T, double p)
{
  if (update(T, p)) return state_->first_partial_deriv(iDmass, iP, iT);
  return 0.0;
};


CoolProp::phases
H2O_DensityCoolProp::get_phase(double T, double p)
{
  if (update(T, p)) return state_->phase();
  return CoolProp::phases::iphase_liquid;
}


bool
H2O_DensityCoolProp::update(double T, double p)
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
