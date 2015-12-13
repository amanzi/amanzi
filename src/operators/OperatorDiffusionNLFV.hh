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

#include <strings.h>

// TPLs
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

/*
  Auxiliary face-based structure for scheme stencils. The stencils 
  are used in the flux calculation. They are not mapped directly 
  on the face-based local matrices that still follow the classical 
  FV structure.
*/

struct FaceStencil {
 public:
  FaceStencil() {};
  ~FaceStencil() {};

  // Use space dimension to allocate enough data for stencils.
  void Init(int d) {
    weights.resize(2 * d);
    stencil.resize(2 * d);
    faces.resize(2 * d);
  }

 public:
  AmanziGeometry::Point p;  // harmonic averaging point
  double gamma;  // coefficient for cell with the lowest id 

  std::vector<double> weights;  // weights in positive decompositions
  std::vector<int> stencil;  // ids of cells in positive decompositions
  std::vector<int> faces;  // ids of interface faces
};

class OperatorDiffusionNLFV : public virtual OperatorDiffusion {
 public:
  OperatorDiffusionNLFV(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusion(global_op),
      stencils_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV;
    InitDiffusion_(plist);
  }

  OperatorDiffusionNLFV(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      stencils_initialized_(false)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV;
    InitDiffusion_(plist);
  }

  OperatorDiffusionNLFV(Teuchos::ParameterList& plist,
                        const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusion(mesh),
      stencils_initialized_(false)
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
  bool stencils_initialized_;
  std::vector<FaceStencil> stencils_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif
