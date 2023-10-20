/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Elasticity operator is used for describing soil deformation or fluid flow (Stokes
and Navier-Stokes).

* `"method`" [string] defines a discretization method. The available
  options are `"BernardiRaugel`".

* `"schema`" [list] defines a discretization schema.

  * `"location`" [Array(string)] defines geometric location of degrees of freedom.

  * `"type`" [Array(string)] defines type of degrees of freedom. The available options
    are `"scalar`" and `"normal component`".

  * `"number`" [Array(int)] indicates how many time this degree of freedom is repeated.

.. code-block:: xml

  <ParameterList name="elasticity operator">
    <Parameter name="method" type="string" value="BernardiRaugel"/>
    <ParameterList name="schema">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_OPERATOR_PDE_ELASTICITY_HH_
#define AMANZI_OPERATOR_PDE_ELASTICITY_HH_

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

// Amanzi::Operators
#include "BilinearForm.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_Elasticity : public PDE_HelperDiscretization {
 public:
  PDE_Elasticity(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh), K_(Teuchos::null), K_default_(1.0)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_ELASTICITY;
    Init_(plist);
  }

  // main virtual members
  // -- setup
  void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor>>& K);
  void SetTensorCoefficient(double K);

  // -- creation of an operator
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;

  // -- modify matrix due to boundary conditions
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial
  //      function, i.e. zeros go in the corresponding matrix columns and
  //      right-hand side is modified using BC values. This is the optional
  //      parameter that enforces symmetry for a symmetric tree operators.
  //    essential_eqn=true indicates that the operator places a positive number on
  //      the main matrix diagonal for the case of essential BCs. This is not
  //      valid for all BCs.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

  // -- postprocessing: calculated stress u from displacement p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override{};

 protected:
  void Init_(Teuchos::ParameterList& plist);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K_;
  double K_default_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;
  AmanziMesh::Entity_kind base_;
};

} // namespace Operators
} // namespace Amanzi

#endif
