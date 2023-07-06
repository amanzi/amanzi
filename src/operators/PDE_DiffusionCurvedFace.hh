/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

// PDE_DiffusionCurvedFace: elliptic operators using the MFD discretization for curved faces


#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_CURVED_FACE_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_CURVED_FACE_HH_

#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"

/*!

Additional options available only for the MFD family of discretizations include:

* `"nonlinear coefficient`" ``[string]`` specifies a method for treating nonlinear
  diffusion coefficient, if any. Available options are `"none`", `"upwind:
  face`", `"divk: cell-face`" (default), `"divk: face`", `"standard: cell`",
  `"divk: cell-face-twin`" and `"divk: cell-grad-face-twin`".  Symmetry
  preserving methods are the divk-family of methods and the classical
  cell-centered method (`"standard: cell`"). The first part of the name
  indicates the base scheme.  The second part (after the semi-column)
  indicates required components of the composite vector that must be provided
  by a physical PK.

* `"discretization secondary`" ``[string]`` specifies the most robust
  discretization method that is used when the primary selection fails to
  satisfy all a priori conditions.  This is typically `"mfd: default`", and is
  used only when an MFD `"discretization primary`" is used.

* `"schema`" ``[Array(string)]`` defines the operator stencil. It is a collection of
  geometric objects.  Typically this is set by the implementation and is not provided.

* `"preconditioner schema`" ``[Array(string)]`` **{face,cell}** Defines the
  preconditioner stencil.  It is needed only when the default assembling
  procedure is not desirable. If skipped, the `"schema`" is used instead.
  In addition to the default, **{face}** may be used, which forms the Schur
  complement.

* `"diffusion tensor`" ``[string]`` allows us to solve problems with symmetric and
  non-symmetric (but positive definite) tensors. Available options are *symmetric*
  (default) and *nonsymmetric*.

*/

/*
  Properties of a minimum norm solution are described here:
  https://faculty.math.illinois.edu/~mlavrov/docs/484-spring-2019/ch4lec4.pdf
*/

namespace Amanzi {
namespace Operators {

class PDE_DiffusionCurvedFace : public virtual PDE_Diffusion {
 public:
  PDE_DiffusionCurvedFace(Teuchos::ParameterList& plist, const Teuchos::RCP<Operator>& global_op)
    : PDE_Diffusion(global_op), plist_(plist), factor_(1.0)
  {
    pde_type_ = PDE_DIFFUSION_MFD_CURVED_FACE;
    Init_(plist);
  }

  PDE_DiffusionCurvedFace(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), plist_(plist), factor_(1.0)
  {
    pde_type_ = PDE_DIFFUSION_MFD_CURVED_FACE;
    Init_(plist);
  }

  PDE_DiffusionCurvedFace(Teuchos::ParameterList& plist, const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), plist_(plist), factor_(1.0)
  {
    pde_type_ = PDE_DIFFUSION_MFD_CURVED_FACE;
    Init_(plist);
  }

  virtual void
  SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K) override;
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) override;

  // -- To calculate elemetal matrices, we can use input parameters flux
  //    and u from the previous nonlinear iteration. Otherwise, use null-pointers.
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  // modify matrix due to boundary conditions
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial
  //      function, i.e. zeros go in the corresponding matrix columns and
  //      right-hand side is modified using BC values. This is the optional
  //      parameter that enforces symmetry for a symmetric tree operators.
  //    essential_eqn=true indicates that the operator places a positive number on
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementation trick.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

  // -- by breaking p-lambda coupling.
  virtual void ModifyMatrices(const CompositeVector& u) override;

  // -- by rescaling mass and stiffness matrices.
  virtual void ScaleMassMatrices(double s) override;
  virtual void ScaleMatricesColumns(const CompositeVector& s) override{};

  // main virtual members after solving the problem
  // -- calculate the flux variable.
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  virtual void UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u,
                                              double scalar_factor = 1.0) override {};

  virtual void UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u,
                                              const Teuchos::Ptr<const CompositeVector>& factor) override {};

  // access
  std::shared_ptr<std::vector<AmanziGeometry::Point>> get_bf() { return bf_; }

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const override { return 0.0; }
  virtual double ComputeGravityFlux(int f) const override { return 0.0; }

  // developer checks
  void set_factor(double factor) { factor_ = factor; }

 protected:
  void Init_(Teuchos::ParameterList& plist);
  void CreateMassMatrices_();
  void LSProblemSetupMatrix_(std::vector<WhetStone::DenseMatrix>& matrices);
  void LSProblemSetupRHS_(CompositeVector& rhs, int i0);
  void LSProblemPrimarySolution_(const CompositeVector& sol, int i0);

 protected:
  Teuchos::ParameterList plist_;
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  bool mass_matrices_initialized_;
  double factor_;

  std::shared_ptr<std::vector<AmanziGeometry::Point>> bf_;

  int schema_prec_dofs_;
};

} // namespace Operators
} // namespace Amanzi


#endif
