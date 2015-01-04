/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for TC implementations.
*/

#include <string>
#include "twophase_thermal_conductivity_factory.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* method for instantiating implementations
****************************************************************** */
Teuchos::RCP<ThermalConductivityTwoPhase> 
   ThermalConductivityTwoPhaseFactory::createThermalConductivityModel(Teuchos::ParameterList& plist)
{
  std::string tc_typename = plist.get<std::string>("thermal conductivity type");
  return Teuchos::rcp(CreateInstance(tc_typename, plist));
};

}  // namespace Energy
}  // namespace Amanzi

