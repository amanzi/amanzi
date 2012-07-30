#ifndef __FLOW_BC_FACTORY_HPP__
#define __FLOW_BC_FACTORY_HPP__

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

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "boundary_function.hh"

namespace Amanzi {
namespace AmanziFlow {

class FlowBCFactory {
 public:
  FlowBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
     : mesh_(mesh), params_(params) {}
  ~FlowBCFactory() {};
  
  BoundaryFunction* createPressure() const;
  BoundaryFunction* createMassFlux() const;
  BoundaryFunction* createStaticHead(double, double, AmanziGeometry::Point&) const;
  BoundaryFunction* createSeepageFace() const;

 private:
  void processPressureList(Teuchos::ParameterList&, BoundaryFunction*) const;
  void processPressureSpec(Teuchos::ParameterList&, BoundaryFunction*) const;
  void processMassFluxList(Teuchos::ParameterList&, BoundaryFunction*) const;
  void processMassFluxSpec(Teuchos::ParameterList&, BoundaryFunction*) const;
  void processStaticHeadList(double p0, double rho, AmanziGeometry::Point& gravity, 
                             Teuchos::ParameterList&, BoundaryFunction*) const;
  void processStaticHeadSpec(double p0, double rho, AmanziGeometry::Point& gravity, 
                             Teuchos::ParameterList&, BoundaryFunction*) const;
  void processSeepageFaceList(Teuchos::ParameterList&, BoundaryFunction*) const;
  void processSeepageFaceSpec(Teuchos::ParameterList&, BoundaryFunction*) const;
     
 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  const Teuchos::RCP<Teuchos::ParameterList>& params_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif // AMANZI_FLOW_BC_FACTORY_HH_
