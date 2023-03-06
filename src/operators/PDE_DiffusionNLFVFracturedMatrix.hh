/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  This implements the PDE_DiffusionNLFV interface for fracture matrix.
  Two DOFS on each fracture face are used.
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_FRACTURED_MATRIX_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_FRACTURED_MATRIX_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

#include "PDE_DiffusionNLFV.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionNLFVFracturedMatrix : public PDE_DiffusionNLFV {
 public:
  PDE_DiffusionNLFVFracturedMatrix(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<Operator>& global_op)
    : PDE_DiffusionNLFV(plist, global_op),
      PDE_Diffusion(global_op),
      stencil_initialized_(false)
  {
    pde_type_ = PDE_DIFFUSION_NLFV_FRACTURED_MATRIX;
    Init_(plist);
  }

  PDE_DiffusionNLFVFracturedMatrix(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_DiffusionNLFV(plist, mesh),
      PDE_Diffusion(mesh),
      stencil_initialized_(false)
  {
    pde_type_ = PDE_DIFFUSION_NLFV_FRACTURED_MATRIX;
    Init_(plist);
  }

  PDE_DiffusionNLFVFracturedMatrix(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
    : PDE_DiffusionNLFV(plist, mesh),
      PDE_Diffusion(mesh),
      stencil_initialized_(false)
  {
    pde_type_ = PDE_DIFFUSION_NLFV_FRACTURED_MATRIX;
    Init_(plist);
  }

  // main virtual members
  // -- setup
  using PDE_Diffusion::Setup;
  virtual void
  SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K) override
  {
    K_ = K;
  }
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) override;

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  virtual void UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u,
                                              double scalar_limiter) override;

  virtual void
  UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                 const Teuchos::Ptr<const CompositeVector>& u,
                                 const Teuchos::Ptr<const CompositeVector>& factor) override;

  // -- after solving the problem: postrocessing
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;
  virtual void ModifyMatrices(const CompositeVector& u) override{};
  virtual void ScaleMassMatrices(double s) override{};
  virtual void ScaleMatricesColumns(const CompositeVector& s) override{};

  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const override { return 0.0; }
  virtual double ComputeGravityFlux(int f) const override { return 0.0; }

 protected:
  // virtual functions for derived clases
  // -- processing of control parameters
  void Init_(Teuchos::ParameterList& plist);
  // -- solution can be modified on boundary faces. This reflects specifics
  //    of nonlinear FV schemes, see implementation in the derived classes.
  // virtual double MapBoundaryValue_(int f, double u) { return u; }

 protected:
  void InitStencils_();
  void OneSidedFluxCorrections_(int i0, const CompositeVector& u, CompositeVector& sideflux);
  void OneSidedWeightFluxes_(int i0, const CompositeVector& u, CompositeVector& sideflux);
  void OneSidedNeumannCorrections_(const CompositeVector& u, CompositeVector& sideflux);
  int OrderCellsByGlobalId_(const AmanziMesh::Entity_ID_List& cells, int& c1, int& c2);
  int NLTPFAContributions_(int f, double& tc1, double& tc2);

 protected:
  int dim_;
  int newton_correction_;

  bool stencil_initialized_;
  Teuchos::RCP<CompositeVector> stencil_data_;
  std::vector<Teuchos::RCP<Epetra_IntVector>> stencil_faces_;
  std::vector<Teuchos::RCP<Epetra_IntVector>> stencil_cells_;
};

} // namespace Operators
} // namespace Amanzi

#endif
