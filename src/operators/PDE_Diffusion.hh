/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

//! Diffusion generates local Ops and global Operators for an elliptic operator.
#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"

/*!

Diffusion is the most frequently used operator. It employs the old schema.

* `"pks operator name`" [list] a PK specific name for the diffusion operator.

  * `"discretization primary`" [string] specifies an advanced discretization method that
    has useful properties under some a priori conditions on the mesh and/or permeability tensor.
    The available options are `"mfd: optimized for sparsity`", `"mfd: optimized for monotonicity`",
    `"mfd: default`", `"mfd: support operator`", `"mfd: two-point flux approximation`",
    `"fv: default`", and `"nlfv: default`".
    The first option is recommended for general meshes.
    The second option is recommended for orthogonal meshes and diagonal absolute 
    permeability tensor. 

  * `"discretization secondary`" [string] specifies the most robust discretization method
    that is used when the primary selection fails to satisfy all a priori conditions.
    Default value is equal to that for the primary discretization.

  * `"diffusion tensor`" [string] specifies additional properties of the diffusion tensor.
    It allows us to solve problems with non-symmetric but positive definite tensors. 
    Available options are *symmetric* (default) and *nonsymmetric*.

  * `"nonlinear coefficient`" [string] specifies a method for treating nonlinear diffusion
    coefficient, if any. Available options are `"none`", `"upwind: face`", `"divk: cell-face`" (default),
    `"divk: face`", `"standard: cell`", and `"divk: cell-face-twin`".
    Symmetry preserving methods are the divk-family of methods and the classical cell-centered
    method (`"standard: cell`"). The first part of the name indicates the base scheme.
    The second part (after the semi-column) indicates required components of the composite vector
    that must be provided by a physical PK.
    Default is `"none`".

  * `"schema`" [Array(string)] defines the operator stencil. It is a collection of 
    geometric objects. It equals to `"{cell}`" for finite volume schemes. 
    It is typically `"{face, cell}`" for mimetic discretizations.

  * `"preconditioner schema`" [Array(string)] defines the preconditioner stencil.
    It is needed only when the default assembling procedure is not desirable. 
    If skipped, the `"schema`" is used instead. 

  * `"gravity`" [bool] specifies if flow is driven also by the gravity.

  * `"gravity term discretization`" [string] selects a model for discretizing the 
    gravity term. Available options are `"hydraulic head`" [default] and `"finite volume`". 
    The first option starts with equation for the shifted solution, i.e. the hydraulic head,
    and derives gravity discretization by the reserve shifting.
    The second option is based on the divergence formula.

  * `"gravity magnitude`" [double] defined magnitude of the gravity vector.

  * `"Newton correction`" [string] specifies a model for correction (non-physical) terms 
    that must be added to the preconditioner. These terms approximate some Jacobian terms.
    Available options are `"true Jacobian`" and `"approximate Jacobian`".
    The FV scheme accepts only the first options. The othre schemes accept only the second option.

  * `"scaled constraint equation`" [bool] rescales flux continuity equations on mesh faces.
    These equations are divided by the nonlinear coefficient. This option allows us to 
    treat the case of zero nonlinear coefficient. At moment this feature does not work 
    with non-zero gravity term. Default is *false*.

  * `"constraint equation scaling cutoff"`" [double] specifies the cutoff value for
    applying rescaling strategy described above.  

  * `"consistent faces`" [list] may contain a `"preconditioner`" and
    `"linear operator`" list (see sections Preconditioners_ and LinearSolvers_
    respectively).  If these lists are provided, and the `"discretization
    primary`" is of type `"mfd: *`", then the diffusion method
    UpdateConsistentFaces() can be used.  This method, given a set of cell
    values, determines the faces constraints that satisfy the constraint
    equation in MFD by assembling and inverting the face-only system.  This is
    not currently used by any Amanzi PKs.

  * `"fracture`" [Array(string)] provides list of regions that defines a fracture network.
    This parameter is used only by the coupled flow PK.


Example:

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
    <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
    <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
    <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
    <Parameter name="gravity" type="bool" value="true"/>
    <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
    <Parameter name="gravity magnitude" type="double" value="9.81"/>
    <Parameter name="nonlinear coefficient" type="string" value="upwind: face"/>
    <Parameter name="Newton correction" type="string" value="true Jacobian"/>

    <ParameterList name="consistent faces">
      <ParameterList name="linear solver">
        ...
      </ParameterList>
      <ParameterList name="preconditioner">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>

This example creates a p-lambda system, i.e. the pressure is
discretized in mesh cells and on mesh faces. 
The preconditioner is defined on faces only, i.e. cell-based unknowns
are eliminated explicitly and the preconditioner is applied to the
Schur complement.

*/


/*
  Ghost elements of composite vectors k_ and dkdp_ are NOT used to be up to date.
  They always have to be updated from master elements before accessing their values.
*/

namespace Amanzi {
namespace Operators {

class PDE_Diffusion : public PDE_HelperDiscretization {
 public:
  PDE_Diffusion(const Teuchos::RCP<Operator>& global_op)
    : PDE_HelperDiscretization(global_op),
      K_(Teuchos::null),
      k_(Teuchos::null),
      dkdp_(Teuchos::null),
      const_k_(1.0){};

  PDE_Diffusion(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh),
      K_(Teuchos::null),
      k_(Teuchos::null),
      dkdp_(Teuchos::null),
      const_k_(1.0){};

  PDE_Diffusion(const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh),
      K_(Teuchos::null),
      k_(Teuchos::null),
      dkdp_(Teuchos::null),
      const_k_(1.0){};

  virtual ~PDE_Diffusion() = default;

  // main virtual members
  // -- setup
  virtual void
  SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K) = 0;
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) = 0;
  void SetConstantScalarCoefficient(double k) { const_k_ = k; }
  void SetConstantTensorCoefficient(const WhetStone::Tensor& K) { const_K_ = K; }

  // -- creation of an operator
  virtual void UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u,
                                              double scalar_factor = 1.0) = 0;

  virtual void
  UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                 const Teuchos::Ptr<const CompositeVector>& u,
                                 const Teuchos::Ptr<const CompositeVector>& factor) = 0;

  // -- matrix modifications
  virtual void ModifyMatrices(const CompositeVector& u) = 0;
  virtual void ScaleMassMatrices(double s) = 0;
  virtual void ScaleMatricesColumns(const CompositeVector& s) = 0;

  // -- default implementation
  virtual void Setup(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K,
                     const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp)
  {
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
  }

  // -- working with consistent faces -- may not be implemented
  virtual int UpdateConsistentFaces(CompositeVector& u)
  {
    Errors::Message msg("Diffusion implementation does not support working with consistent faces.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // additional interface
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const = 0;
  virtual double ComputeGravityFlux(int f) const = 0;

  // access
  int schema_prec_dofs() { return global_op_schema_; }
  int schema_dofs() { return local_op_schema_; }

  Teuchos::RCP<const Op> jacobian_op() const { return jac_op_; }
  Teuchos::RCP<Op> jacobian_op() { return jac_op_; }
  void set_jacobian_op(const Teuchos::RCP<Op>& op);
  int schema_jacobian() { return jac_op_schema_; }

  int little_k() const { return little_k_; }
  CompositeVectorSpace little_k_space() const
  {
    CompositeVectorSpace out;
    out.SetMesh(mesh_);
    out.SetGhosted();
    if (little_k_ == OPERATOR_LITTLE_K_NONE) { return out; }
    if (little_k_ != OPERATOR_LITTLE_K_UPWIND) { out.AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1); }
    if (little_k_ != OPERATOR_LITTLE_K_STANDARD) { out.AddComponent("face", AmanziMesh::Entity_kind::FACE, 1); }
    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN) { out.AddComponent("twin", AmanziMesh::Entity_kind::FACE, 1); }
    return out;
  }

 protected:
  // -- additional interface on non-manifolds
  virtual void UpdateFluxManifold_(const Teuchos::Ptr<const CompositeVector>& u,
                                   const Teuchos::Ptr<CompositeVector>& flux);

 protected:
  Teuchos::RCP<const std::vector<WhetStone::Tensor>> K_;
  WhetStone::Tensor const_K_;
  bool K_symmetric_;

  // nonlinear coefficient and its representation
  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  double const_k_;
  int little_k_;

  // additional operators
  Teuchos::RCP<Op> jac_op_;
  int global_op_schema_, local_op_schema_, jac_op_schema_;
};


/* ******************************************************************
* Default implementations
****************************************************************** */
inline void 
PDE_Diffusion::UpdateFluxManifold_(const Teuchos::Ptr<const CompositeVector>& u,
                                   const Teuchos::Ptr<CompositeVector>& flux)
{
  Errors::Message msg("Diffusion does not know how to compute flux on manifolds.");
  Exceptions::amanzi_throw(msg);
}


/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
inline void
PDE_Diffusion::set_jacobian_op(const Teuchos::RCP<Op>& op)
{
  if (global_operator().get()) {
    if (local_op().get()) {
      auto index = std::find(global_operator()->begin(), global_operator()->end(), jac_op_) -
                   global_operator()->begin();
      if (index != global_operator()->size())
        global_operator()->OpPushBack(op);
      else
        global_operator()->OpReplace(op, index);
    } else {
      global_operator()->OpPushBack(op);
    }
  }
  jac_op_ = op;
}

} // namespace Operators
} // namespace Amanzi

#endif
