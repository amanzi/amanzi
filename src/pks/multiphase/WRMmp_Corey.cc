/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

#include "MultiphaseDefs.hh"
#include "WRMmp_Corey.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRMmp_Corey::WRMmp_Corey(Teuchos::ParameterList& plist)
{
  srl_ = plist.get<double>("residual saturation liquid", 0.0);
  srg_ = plist.get<double>("residual saturation gas", 0.0);
  pcap_ = plist.get<double>("capillary pressure");
}


/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRMmp_Corey::k_relative(double sl, int phase)
{
  double sle = (sl - srl_ - srg_) / (1.0 - srl_ - srg_);
  if (phase == MULTIPHASE_PHASE_LIQUID) {
    return std::pow(sle, 4.0);
  }
  else if (phase == MULTIPHASE_PHASE_GAS) {
    return std::pow(1.0 - sle, 2.0) * (1.0 - sle * sle);
  }

  return 0.0;
}


/* ******************************************************************
* Derivative of relative permeability wrt liquid saturation. 
****************************************************************** */
double WRMmp_Corey::dKdS(double sl, int phase)
{
  double factor = 1.0 / (1.0 - srl_ - srg_);
  double sle = (sl - srl_) / (1.0 - srl_ - srg_);
  if (phase == MULTIPHASE_PHASE_LIQUID) {
    return 4.0 * factor * std::pow(sle, 3.0);
  }
  else if (phase == MULTIPHASE_PHASE_GAS) {
    return -2.0 * factor * std::pow(1.0 - sle, 2.0) * (1.0 + 2.0 * sle);
  }
}


/* ******************************************************************
* Capillary pressure formula.
****************************************************************** */
double WRMmp_Corey::capillaryPressure(double sl)
{
  if (sl <= srl_) return pcap_;
  return pcap_ * (1.0 - sl) / (1.0 - srl_);
}


/* ******************************************************************
* Derivative of capillary pressure. Hard-coded Brooks-Corey
* with Pd = 1, gamma = 3. Assume the saturation is of the wetting phase
****************************************************************** */
double WRMmp_Corey::dPc_dS(double sl)
{
  if (sl <= srl_) return 0.0;
  return -pcap_ / (1.0 - srl_);
}

}  // namespace Multiphase
}  // namespace Amanzi

