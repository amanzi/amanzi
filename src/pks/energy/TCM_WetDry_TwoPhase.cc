/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  Simple model of two-phase thermal conductivity, based upon:
    - Interpolation between saturated and dry conductivities via a Kersten number.
    - Power-law Kersten number.
  See ATS process model documentation's permafrost model for details.
*/

#include <cmath>
#include "TCM_WetDry_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Simple constructor
****************************************************************** */
TCM_WetDry_TwoPhase::TCM_WetDry_TwoPhase(Teuchos::ParameterList& plist) : plist_(plist)
{
  InitializeFromPlist_();
}


/* ******************************************************************
* Equation of state
****************************************************************** */
double
TCM_WetDry_TwoPhase::ThermalConductivity(double poro, double sat_liq)
{
  double kersten = pow(sat_liq + eps_, alpha_);
  return k_dry_ + (k_wet_ - k_dry_) * kersten;
};


/* ******************************************************************
* Initialization
****************************************************************** */
void
TCM_WetDry_TwoPhase::InitializeFromPlist_()
{
  eps_ = plist_.get<double>("epsilon", 1.e-10);
  alpha_ = plist_.get<double>("unsaturated alpha", 1.0);
  k_dry_ = plist_.get<double>("thermal conductivity, dry");
  k_wet_ = plist_.get<double>("thermal conductivity, wet");
};

} // namespace Energy
} // namespace Amanzi
