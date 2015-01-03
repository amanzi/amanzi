/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for EOS implementations.
*/

#ifndef PK_ENERGY_RELATIONS_TC_TWOPHASE_FACTORY_HH_
#define PK_ENERGY_RELATIONS_TC_TWOPHASE_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "twophase_thermal_conductivity.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityTwoPhaseFactory : public Utils::Factory<ThermalConductivityTwoPhase> {
 public:
  Teuchos::RCP<ThermalConductivityTwoPhase> createThermalConductivityModel(Teuchos::ParameterList& plist);
};

}  // namespace Energy
}  // namespace Amanzi

#endif
