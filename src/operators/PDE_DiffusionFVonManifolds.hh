/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  The FV scheme for diffusion equation on manifolds. The key difference
  with another FV implemnetation is the flux continuity condition.
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FV_MANIFOLDS_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FV_MANIFOLDS_HH_

#include <strings.h>

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Opertors
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionFVonManifolds : public virtual PDE_Diffusion, public PDE_DiffusionWithGravity {
 public:
  PDE_DiffusionFVonManifolds(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<Operator>& global_op,
                             bool flag)
    : PDE_Diffusion(global_op), PDE_DiffusionWithGravity(global_op), beta_initialized_(false)
  {
    pde_type_ = PDE_DIFFUSION_FV_MANIFOLDS;
    Init_(plist);
  }

  PDE_DiffusionFVonManifolds(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), PDE_DiffusionWithGravity(mesh), beta_initialized_(false)
  {
    pde_type_ = PDE_DIFFUSION_FV_MANIFOLDS;
    Init_(plist);
  }

  // main virtual members
  // -- setup
  using PDE_Diffusion::Setup;
  virtual void
  SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K) override;
  // -- only one component is used:
  //    face component if present - use it, ignore other components
  //    only cell component is present - use it.
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) override;

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;
  virtual void UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                              const Teuchos::Ptr<const CompositeVector>& u,
                                              double scalar_factor = 1.0) override{};

  virtual void
  UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                 const Teuchos::Ptr<const CompositeVector>& u,
                                 const Teuchos::Ptr<const CompositeVector>& factor) override{};

  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // -- modify the local operator
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;
  virtual void ModifyMatrices(const CompositeVector& u) override{};
  virtual void ScaleMassMatrices(double s) override{};
  virtual void ScaleMatricesColumns(const CompositeVector& s) override{};

  // -- transmisibility is multi-valued for manifolds and is skipped
  virtual double ComputeTransmissibility(int f) const override { return 0.0; }
  virtual double ComputeGravityFlux(int f) const override { return 0.0; }

 protected:
  void ComputeBeta_();

  void Init_(Teuchos::ParameterList& plist);

 private:
  Teuchos::RCP<CompositeVector> beta_; // static part of transmissibility
  bool beta_initialized_, gravity_;
};

} // namespace Operators
} // namespace Amanzi


#endif
