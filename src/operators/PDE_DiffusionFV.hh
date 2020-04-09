/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Daniil Svyatskiy (dasvyat@lanl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FV_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FV_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Opertors
#include "PDE_Diffusion.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionFV : public virtual PDE_Diffusion {
 public:
  PDE_DiffusionFV(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<Operator>& global_op)
    : PDE_Diffusion(global_op), transmissibility_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV;
    Init_(plist);
  }

  PDE_DiffusionFV(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), transmissibility_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV;
    Init_(plist);
  }

  PDE_DiffusionFV(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), transmissibility_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_FV;
    Init_(plist);
  }

  // main virtual members
  // -- setup
  using PDE_Diffusion::Setup;
  virtual void SetTensorCoefficient(
    const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K) override;
  virtual void SetScalarCoefficient(
    const Teuchos::RCP<const CompositeVector>& k,
    const Teuchos::RCP<const CompositeVector>& dkdp) override;

  // -- create an operator
  virtual void
  UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                 const Teuchos::Ptr<const CompositeVector>& u) override;
  virtual void UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    double scalar_factor = 1.0) override;

  virtual void UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor) override;

  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;
  virtual void
  UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
                        const Teuchos::Ptr<CompositeVector>& flux) override;

  // -- modify an operator
  virtual void
  ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;
  virtual void ModifyMatrices(const CompositeVector& u) override{};
  virtual void ScaleMassMatrices(double s) override{};

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const override;
  virtual double ComputeGravityFlux(int f) const override { return 0.0; }

  // access
  const CompositeVector& transmissibility() { return *transmissibility_; }

 protected:
  void ComputeTransmissibility_();

  void AnalyticJacobian_(const CompositeVector& solution);

  virtual void
  ComputeJacobianLocal_(int mcells, int f, int face_dir_0to1, int bc_model,
                        double bc_value, double* pres, double* dkdp_cell,
                        WhetStone::DenseMatrix& Jpp);

  void Init_(Teuchos::ParameterList& plist);

 protected:
  Teuchos::RCP<CompositeVector> transmissibility_;
  bool transmissibility_initialized_;

  int newton_correction_;
  bool exclude_primary_terms_;
};

} // namespace Operators
} // namespace Amanzi


#endif
