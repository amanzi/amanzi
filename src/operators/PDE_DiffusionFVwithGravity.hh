/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_FV_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_FV_WITH_GRAVITY_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Operators
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionFVwithGravity : public PDE_DiffusionFV {
 public:
  PDE_DiffusionFVwithGravity(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<Operator>& global_op) :
      PDE_DiffusionFV(plist, global_op)
  {}

  PDE_DiffusionFVwithGravity(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_DiffusionFV(plist, mesh)
  {}
  
  virtual void Init() override;

  // main virtual members
  // -- setup
  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho) override;
  virtual void SetDensity(double rho) override;

  // -- create a lineratized operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;
  virtual void UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
                                     const Teuchos::Ptr<CompositeVector>& flux) override;

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;
  virtual void ModifyMatrices(const CompositeVector& u) override {};
  virtual void ScaleMassMatrices(double s) override {
    ComputeTransmissibility_(gravity_term_);
    transmissibility_->Scale(s);
  };

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const override;

  // access
  const CompositeVector& gravity_terms() { return *gravity_term_; }

 protected:
  virtual void ComputeJacobianLocal_(
      int mcells, int f, int face_dir_0to1, int bc_model, double bc_value,
      double *pres, double *dkdp_cell, WhetStone::DenseMatrix& Jpp) override;

  void ComputeTransmissibility_(Teuchos::RCP<CompositeVector> g_cv);

 protected:
  Teuchos::RCP<CompositeVector> gravity_term_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
