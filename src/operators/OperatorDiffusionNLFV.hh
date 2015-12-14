/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  OperatorDiffusionNLFV implements the OperatorDiffusion interface 
  using nonlinear finite volumes.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_NLFV_HH_
#define AMANZI_OPERATOR_DIFFUSION_NLFV_HH_

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
#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionNLFV : public virtual OperatorDiffusion {
 public:
  OperatorDiffusionNLFV(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusion(global_op),
      stencil_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV;
    InitDiffusion_(plist);
  }

  OperatorDiffusionNLFV(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      stencil_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV;
    InitDiffusion_(plist);
  }

  OperatorDiffusionNLFV(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      stencil_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV;
    InitDiffusion_(plist);
  }

  // main virtual members
  // -- setup
  using OperatorDiffusion::Setup;
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) { K_ = K; }
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) {
    k_ = k;
    dkdp_ = dkdp;
  }

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux) {};

  // -- modify an operator
  virtual void ApplyBCs(bool primary, bool eliminate);
  virtual void ModifyMatrices(const CompositeVector& u) {};
  virtual void ScaleMassMatrices(double s) {};

  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const { return 0.0; }
  virtual double ComputeGravityFlux(int f) const { return 0.0; }

 private:
  virtual void InitDiffusion_(Teuchos::ParameterList& plist);
  void InitStencils_();
  double OneSidedFluxCorrection_(int c, int f, const Epetra_MultiVector& uc, int k);
  
 private:
  int dim_;
  bool stencil_initialized_;
  Teuchos::RCP<CompositeVector> stencil_data_;
  std::vector<Teuchos::RCP<Epetra_IntVector> > stencil_faces_;
  std::vector<Teuchos::RCP<Epetra_IntVector> > stencil_cells_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
