/*
  Energy
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_ENERGY_BC_FACTORY_HH_
#define AMANZI_ENERGY_BC_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "VerboseObject.hh"

#include "EnergyBoundaryFunction.hh"

namespace Amanzi {
namespace Energy {

class EnergyBCFactory {
 public:
  EnergyBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                const Teuchos::RCP<Teuchos::ParameterList>& plist);
  ~EnergyBCFactory() {}
  
  EnergyBoundaryFunction* CreateTemperature(std::vector<int>& submodel) const;
  EnergyBoundaryFunction* CreateEnergyFlux(std::vector<int>& submodel) const;

 private:
  void ProcessTemperatureList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const;
  void ProcessTemperatureSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const;

  void ProcessEnergyFluxList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const;
  void ProcessEnergyFluxSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, EnergyBoundaryFunction* bc) const;

  void PopulateSubmodelFlag(
      const std::vector<std::string>& regions, int flag, std::vector<int>& submodel) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  const Teuchos::RCP<Teuchos::ParameterList>& plist_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
