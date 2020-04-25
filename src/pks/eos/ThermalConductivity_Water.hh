/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Thermal conductivity for liquid water.
*/

#ifndef AMANZI_EOS_THERMAL_CONDUCTIVITY_WATER_HH_
#define AMANZI_EOS_THERMAL_CONDUCTIVITY_WATER_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class ThermalConductivity_Water {
 public:
  explicit
  ThermalConductivity_Water(Teuchos::ParameterList& eos_plist);

  virtual double ThermalConductivity(double T);
  virtual double DThermalConductivityDT(double T);

 protected:
  virtual void InitializeFromPlist_();

 protected:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded
  double ka0_, ka1_, ka2_;
  double kref_, Tref_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
