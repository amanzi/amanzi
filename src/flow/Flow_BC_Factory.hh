/*
This is the flow component of the Amanzi code.

 
Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __FLOW_BC_FACTORY_HH__
#define __FLOW_BC_FACTORY_HH__

/* transient-flow-bc-list(NAME) is:
  <ParameterList name=NAME>
    transient-flow-condition-list
    ...
    transient-flow-condition-list
  </ParameterList>

where transient-flow-condition-list is one of:
    transient-flow-pressure-list
    transient-flow-mass-flux-list
    transient-flow-static-head-list
  
The parameter list name string NAME is arbitrary, and meaningful only 
to the parent parameter list.
Each sublist defines one type of condition for a transient flow problem.
Each specific type of sublist can appear at most once.
This parameter list is given to a boundary condition "factory" which has
methods for instantiating the appropriate boundary condition objects.

1. transient-flow-pressure-list is:
  <ParameterList name="pressure">
    transient-flow-pressure-spec(NAME_1)
    ...
    transient-flow-pressure-spec(NAME_N)
  </ParameterList>
  
2. transient-flow-mass-flux-list is:
  <ParameterList name="mass flux">
    transient-flow-mass-flux-spec(NAME_1)
    ...
    transient-flow-mass-flux-spec(NAME_N)
  </ParameterList>
  
3. transient-flow-static-head-list is:
  <ParameterList name="static head">
    transient-flow-static-head-spec(NAME_1)
    ...
    transient-flow-static-head-spec(NAME_N)
  </ParameterList>
  
Each spec sublist defines one part of the of total boundary condition of that type.
The name strings NAME_1, ..., NAME_N are arbitrary but must be unique
within the parent parameter list.  They may be used in error messages
or diagnostic logging.  Example: "north cribs"
  Each of the following kinds of spec parameter lists have a "regions"
parameter that is an array of 1 or more region names that specify the
portion of the mesh boundary where the condition is to be applied.
    
transient-flow-pressure-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="regions" type="Array string" value=string-array />
    function-factory-list("boundary pressure")
  </ParameterList>
  
The function-factory-list should define a function whose argument
is the vector (t, x, y, z).
    
transient-flow-mass-flux-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="Regions" type="Array string" value=string-array />
    function-factory-list("outward mass flux")
  </ParameterList>
  
The function-factory-list should define a function whose argument
is the vector (t, x, y, z).

transient-flow-static-head-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="Regions" type="Array string" value=string-array />
    function-factory-list("water table elevation")
  </ParameterList>
  
The function-factory-list should define a function h whose argument
is the vector (t, x, y); the water table at time t is taken to be
the surface (x, y, h(t,x,y)).
The gravitational acceleration is assumed to be directed in the negative z-direction.
*/

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "BoundaryFunction.hh"


namespace Amanzi {
namespace AmanziFlow {

class FlowBCFactory {
 public:
  FlowBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
     : mesh_(mesh), params_(params) {}
  ~FlowBCFactory() {};
  
  BoundaryFunction* CreatePressure(std::vector<int>& submodel) const;
  BoundaryFunction* CreateMassFlux(std::vector<int>& submodel) const;
  BoundaryFunction* CreateStaticHead(
      double p0, double rho, const AmanziGeometry::Point& gravity, std::vector<int>& submodel) const;
  BoundaryFunction* CreateSeepageFace(std::vector<int>& submodel) const;

 private:
  void ProcessPressureList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;
  void ProcessPressureSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;

  void ProcessMassFluxList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;
  void ProcessMassFluxSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;

  void ProcessSeepageFaceList(
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;
  void ProcessSeepageFaceSpec(
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;

  void ProcessStaticHeadList(
      double p0, double rho, const AmanziGeometry::Point& gravity, 
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;
  void ProcessStaticHeadSpec(
      double p0, double rho, const AmanziGeometry::Point& gravity,
      Teuchos::ParameterList& list, std::vector<int>& submodel, BoundaryFunction* bc) const;
     
  void PopulateSubmodelFlag(
      const std::vector<std::string>& regions, int flag, std::vector<int>& submodel) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  const Teuchos::RCP<Teuchos::ParameterList>& params_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
