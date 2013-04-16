#ifndef TESTS_BC_FACTORY_HH_
#define TESTS_BC_FACTORY_HH_

/* -------------------------------------------------------------------------
ATS

Author: ...
    Ethan Coon (ATS version) (ecoon@lanl.gov)

   transient-flow-bc-list(NAME) is:
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
namespace TestPKs {

class BCFactory {

public:
  BCFactory(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                const Teuchos::ParameterList& plist)
     : mesh_(mesh), plist_(plist) {}

  Teuchos::RCP<Functions::BoundaryFunction> CreateDirichlet() const;
  Teuchos::RCP<Functions::BoundaryFunction> CreateNeumann() const;

private:
  void ProcessDirichletList(const Teuchos::ParameterList&,
                           const Teuchos::RCP<Functions::BoundaryFunction>&) const;
  void ProcessDirichletSpec(const Teuchos::ParameterList&,
                           const Teuchos::RCP<Functions::BoundaryFunction>&) const;
  void ProcessNeumannList(const Teuchos::ParameterList&,
                           const Teuchos::RCP<Functions::BoundaryFunction>&) const;
  void ProcessNeumannSpec(const Teuchos::ParameterList&,
                           const Teuchos::RCP<Functions::BoundaryFunction>&) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  Teuchos::ParameterList plist_;
};

}  // namespace
}  // namespace

#endif // AMANZI_FLOW_BC_FACTORY_HH_
