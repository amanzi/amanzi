/*
  This is the multiphase component of the Amanzi code.
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MULTIPHASE_BC_FACTORY_HH_
#define AMANZI_MULTIPHASE_BC_FACTORY_HH_

/* ******************************************************************
transient-flow-bc-list(NAME) is:
  <ParameterList name=NAME>
    transient-flow-condition-list
    ...
    transient-flow-condition-list
  </ParameterList>

where transient-flow-condition-list is one of:
    transient-flow-saturation-list
    transient-flow-mass-flux-list

Currently, mass flux is bc factory is pulled from flow pk.
  
The parameter list name string NAME is arbitrary, and meaningful only 
to the parent parameter list.
Each sublist defines one type of condition for a transient flow problem.
Each specific type of sublist can appear at most once.
This parameter list is given to a boundary condition "factory" which has
methods for instantiating the appropriate boundary condition objects.

1. transient-flow-saturation-list is:
  <ParameterList name="saturation">
    transient-flow-saturation-spec(NAME_1)
    ...
    transient-flow-saturation-spec(NAME_N)
  </ParameterList>
  
2. transient-flow-mass-flux-list is:
  <ParameterList name="mass flux">
    transient-flow-mass-flux-spec(NAME_1)
    ...
    transient-flow-mass-flux-spec(NAME_N)
  </ParameterList>
  
  
Each spec sublist defines one part of the of total boundary condition of that type.
The name strings NAME_1, ..., NAME_N are arbitrary but must be unique
within the parent parameter list.  They may be used in error messages
or diagnostic logging.  Example: "north cribs"
  Each of the following kinds of spec parameter lists have a "regions"
parameter that is an array of 1 or more region names that specify the
portion of the mesh boundary where the condition is to be applied.
    
transient-flow-saturation-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="regions" type="Array string" value=string-array />
    function-factory-list("boundary saturation")
  </ParameterList>
  
The function-factory-list should define a function whose argument
is the vector (t, x, y, z).
    
transient-flow-mass-flux-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="Regions" type="Array string" value=string-array />
    function-factory-list("outward mass flux")
  </ParameterList>
  
* *******************************************************************/

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "VerboseObject.hh"

#include "FlowBoundaryFunction.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseBCFactory {
 public:
  MultiphaseBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                const Teuchos::RCP<Teuchos::ParameterList>& plist);
  ~MultiphaseBCFactory() { delete vo_; }
  
  Flow::FlowBoundaryFunction* CreateSaturation(std::vector<int>& submodel) const;
  Flow::FlowBoundaryFunction* CreatePressure(std::vector<int>& submodel) const;
  Flow::FlowBoundaryFunction* CreateMolarFraction(std::vector<int>& submodel, int phase) const; 
  Flow::FlowBoundaryFunction* CreateMassFlux(std::vector<int>& submodel, int phase) const;
  Flow::FlowBoundaryFunction* CreateMassFlux(std::vector<int>& submodel) const;
  Flow::FlowBoundaryFunction* CreateHydrogenDensity(std::vector<int>& submodel) const;  

 private:
  void ProcessHydrogenDensityList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessHydrogenDensitySpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessSaturationList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessSaturationSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessPressureList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessPressureSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
   void ProcessMolarFractionList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const; 
  void ProcessMolarFractionSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessMassFluxList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void ProcessMassFluxSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, Flow::FlowBoundaryFunction* bc) const;
  void PopulateSubmodelFlag(
      const std::vector<std::string>& regions, int flag, std::vector<int>& submodel) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  const Teuchos::RCP<Teuchos::ParameterList>& plist_;

  VerboseObject* vo_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
