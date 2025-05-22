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

The diffusion operator supports multiple discretization methods to implement
diffusion.

.. _pde-diffusion-spec:
.. admonition:: pde-diffusion-spec

   * `"discretization primary`" ``[string]`` specifies an advanced discretization method that
     has useful properties under some a priori conditions on the mesh and/or permeability tensor.
     The available options are `"mfd: optimized for sparsity`", `"mfd: optimized for monotonicity`",
     `"mfd: default`", `"mfd: support operator`", `"mfd: two-point flux approximation`",
     `"fv: default`", and `"nlfv: default`".
     The first option is recommended for general meshes.
     The second option is recommended for orthogonal meshes and diagonal absolute
     permeability tensor.

   * `"discretization secondary`" ``[string]`` **optional** specifies the most
     robust discretization method that is used when the primary selection fails
     to satisfy all a priori conditions.  Default value is equal to that for
     the primary discretization.

   * `"diffusion tensor`" ``[string]`` **optional** specifies additional
     properties of the diffusion tensor.  It allows us to solve problems with
     non-symmetric but positive definite tensors.  Available options are
     `"symmetric`" and `"nonsymmetric`".

   * `"schema`" ``[Array(string)]`` **optional** defines the operator
     stencil. It is a collection of geometric objects. It is `"{cell}`" for
     finite volume schemes.  It is typically `"{face, cell}`" for mimetic
     discretizations.  This is typically not set by the user.

   * `"preconditioner schema`" ``[Array(string)]`` **optional** defines the
     preconditioner stencil.  It is needed only when the default assembling
     procedure is not desirable.  If skipped, the `"schema`" is used instead.
     For instance, if `"{face}`" is provided for a `"{face,cell}`" schema, then
     the Schur complement is implemented for the preconditioner.

   * `"gravity`" ``[bool]`` **optional** specifies if flow is driven also by
     the gravity.  This is typically set by the PK using this operator.

   * `"Newton correction`" ``[string]`` **none** specifies a model for
     correction terms that may be added to the preconditioner for nonlinear
     coefficients. These terms approximate some Jacobian terms.  Available
     options are `"true Jacobian`" and `"approximate Jacobian`".  The FV scheme
     accepts only the first option. The other schemes accept only the second
     option.

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
  virtual void ScaleMatricesColumns(const CompositeVector& s);

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
  CompositeVectorSpace little_k_space() const;

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
* Scale face-based matrices.
****************************************************************** */
inline void
PDE_Diffusion::ScaleMatricesColumns(const CompositeVector& s)
{
  if (!(local_op_schema_ & OPERATOR_SCHEMA_BASE_FACE)) AMANZI_ASSERT(false);
  if (!s.HasComponent("cell")) AMANZI_ASSERT(false);

  const auto& s_c = *s.ViewComponent("cell");

  for (int f = 0; f < nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];

    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    for (int n = 0; n < ncells; ++n) {
      double factor = s_c[0][cells[n]];
      for (int m = 0; m < ncells; ++m) { Aface(m, n) *= factor; }
    }
  }
}


/* ******************************************************************
* Default implementations
****************************************************************** */
inline CompositeVectorSpace
PDE_Diffusion::little_k_space() const
{
  CompositeVectorSpace out;
  out.SetMesh(mesh_);
  out.SetGhosted();
  if (little_k_ == OPERATOR_LITTLE_K_NONE) { return out; }
  if (little_k_ != OPERATOR_LITTLE_K_UPWIND) {
    out.AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }
  if (little_k_ != OPERATOR_LITTLE_K_STANDARD) {
    out.AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  }
  if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN) {
    out.AddComponent("twin", AmanziMesh::Entity_kind::FACE, 1);
  }
  return out;
}


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
