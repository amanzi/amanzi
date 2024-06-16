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

.. admonition:: elasticity_op-spec

  * `"method`" ``[string]`` defines a discretization method. The available
    options are `"BernardiRaugel`".

  * `"schema`" ``[list]`` defines a discretization schema.

    * `"location`" ``[Array(string)]`` defines geometric location of degrees of freedom.

    * `"type`" ``[Array(string)]`` defines type of degrees of freedom. The available options
      are `"scalar`" and `"normal component`".

    * `"number`" ``[Array(int)]`` indicates how many time this degree of freedom is repeated.

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
    : PDE_HelperDiscretization(mesh), K_(Teuchos::null)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_ELASTICITY;
  }

  PDE_Elasticity(Teuchos::ParameterList& plist, const Teuchos::RCP<Operator>& global_op)
    : PDE_HelperDiscretization(global_op)
  {
    pde_type_ = PDE_ELASTICITY;
  }

  // setup
  void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor>>& C);
  void SetTensorCoefficient(const WhetStone::Tensor& C);
  void SetScalarCoefficient(const CompositeVector& C);

  // --- Young modulus and Poisson ratio
  void SetTensorCoefficientEnu(const Teuchos::RCP<const CompositeVector>& E,
                               const Teuchos::RCP<const CompositeVector>& nu);
  // --- Shear modulus and bulk modulus
  void SetTensorCoefficientGK(const Teuchos::RCP<const CompositeVector>& G,
                              const Teuchos::RCP<const CompositeVector>& K);

  // main virtual members
  virtual void Init(Teuchos::ParameterList& plist);

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

  // -- postprocessing: place holder
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override{};

  // -- processing: calculate stesses from displacement
  void ComputeHydrostaticStress(const CompositeVector& u, CompositeVector& p);
  void ComputeVolumetricStrain(const CompositeVector& u, CompositeVector& e);

  // -- cell-based algorithms
  WhetStone::Tensor ComputeCellStrain(const CompositeVector& u, int c);

 protected:
  WhetStone::Tensor computeElasticityTensorEnu_(int c);
  WhetStone::Tensor computeElasticityTensorGK_(int c);

 private:
  void ApplyBCs_Kinematic_(const BCs& bc, bool primary, bool eliminate, bool essential_eqn);
  void ApplyBCs_ShearStress_(const BCs& bc, bool primary, bool eliminate, bool essential_eqn);
  void ApplyBCs_Traction_(const BCs& bc, bool primary, bool eliminate, bool essential_eqn);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor>> C_;
  WhetStone::Tensor C_default_;
  Teuchos::RCP<const CompositeVector> E_, nu_;
  Teuchos::RCP<const CompositeVector> G_, K_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;
  AmanziMesh::Entity_kind base_;
};

} // namespace Operators
} // namespace Amanzi

#endif
