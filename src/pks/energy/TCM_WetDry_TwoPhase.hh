/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Simple model of two-phase thermal conductivity, based upon:
   - Interpolation between saturated and dry conductivities via a Kersten number.
   - Power-law Kersten number.
  See ATS process model documentation's permafrost model for details.

  Usage:

    <ParameterList name="thermal conductivity model">
      <Parameter name="thermal conductivity type" type="string" value="two-phase wet/dry"/>

      <Parameter name="thermal conductivity, wet" type="double" value=""/>
      <Parameter name="thermal conductivity, dry" type="double" value=""/>

      <Parameter name="epsilon" type="double" value="1.e-10"/>
      <Parameter name="unsaturated alpha" type="double" value="1.0"/>
    </ParameterList>

  Units: ????
*/

#ifndef PK_ENERGY_TCM_WETDRY_TWOPHASE_HH_
#define PK_ENERGY_TCM_WETDRY_TWOPHASE_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "TCM_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

class TCM_WetDry_TwoPhase : public TCM_TwoPhase {
 public:
  TCM_WetDry_TwoPhase(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_wet_;
  double k_dry_;

 private:
  static Utils::RegisteredFactory<TCM_TwoPhase,TCM_WetDry_TwoPhase> factory_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
