/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*!

A high-order advection operator may have different domain and range and therefore requires two schemas.
The structure of the new schema is described in the previous section.
A high-order advection operator has two terms in a weak formulation, corresponding to 
volume and surface integrals. These two terms are discretixed using two operators with
matrix of types *advection* and *flux*, respectively.


* `"pks operator name`" [list] a PK specific name for the advection operator.

  * `"method`" [string] defines a discretization method. The available option is `"dg modal`".

  * `"method order`" [int] defines method order. For example, the classical low-order finite 
    volume scheme is equivalent to DG of order 0.

  * `"matrix type`" [string] defines matrix type. The supported options are `"advection`"
    and `"flux`".

  * `"dg basis`" [string] defines bases for DG schemes. The available options are 
    `"regularized`" (recommended), `"normalized`", `"orthonormalized`", and `"natural`" 
    (not recommended).

  * `"gradient operator on test function`" [bool] defines place of the gradient operator.
    For integration by parts schemes, the gradient is transfered to a test function.
    This option is needed for discretizing volumetric integrals.

  * `"jump operator on test function`" [bool] defines place of the jump operator.
    For integration by parts schemes, the jump operator is applied to a test function.
    This option is needed for discretizing surface fluxes.

  * `"flux formula`" [string] defines type of the flux. The available options 
    are `"Rusanov`" (default), `"upwind`", `"downwind`", and `"NavierStokes`".

  * `"schema domain`" [list] defines a discretization schema for the operator domain.

  * `"schema range`" [list] defines a discretization schema for the operator range. 

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="method" type="string" value="dg modal"/>
    <Parameter name="method order" type="int" value="2"/>
    <Parameter name="flux formula" type="string" value="Rusanov"/>
    <Parameter name="matrix type" type="string" value="flux"/>
    <Parameter name="jump operator on test function" type="bool" value="true"/>

    <ParameterList name="schema domain">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
    <ParameterList name="schema range">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{cell}"/>
      <Parameter name="type" type="Array(string)" value="{scalar}"/>
      <Parameter name="number" type="Array(int)" value="{1}"/>
    </ParameterList>
  </ParameterList>

In this example, we construct an operator for volumetric integrals in a weak formulation
of advection problem.

The only low-order advection operator in Amanzi is the upwind operator. 
It employes the old schema.

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="base" type="string" value="face"/>
    <Parameter name="schema" type="Array(string)" value="{cell}"/>
    <Parameter name="method order" type="int" value="0"/>
    <Parameter name="matrix type" type="string" value="advection"/>
  </ParameterList>

*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_HH_

#include "PDE_HelperDiscretization.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"


namespace Amanzi {
namespace Operators {

class PDE_Advection : public PDE_HelperDiscretization {
 public:
  PDE_Advection(Teuchos::ParameterList& plist, const Teuchos::RCP<Operator>& global_op)
    : PDE_HelperDiscretization(global_op){};

  PDE_Advection(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh)
  {
    global_op_ = Teuchos::null;
  }

  virtual ~PDE_Advection(){};

  // main members
  // -- setup
  virtual void Setup(const CompositeVector& u) = 0;

  // -- standard interface for flux calculation
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override{};

  // -- extended interface for flux calculation
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h,
                          const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::RCP<BCs>& bc,
                          const Teuchos::Ptr<CompositeVector>& u) = 0;

 protected:
  std::string name_;
};

} // namespace Operators
} // namespace Amanzi

#endif
