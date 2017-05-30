/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete gravity operator blended with the MFD diffusion operator.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_MFD_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_MFD_WITH_GRAVITY_HH_

#include "Epetra_IntVector.h"

#include "Tensor.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

#include "DiffusionMFD.hh"
#include "DiffusionWithGravity.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class DiffusionMFDwithGravity : public DiffusionMFD,
                                public DiffusionWithGravity {
 public:
  DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<Operator>& global_op) :
      DiffusionMFD(plist, global_op),
      DiffusionWithGravity(global_op),
      Diffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
  }
  
  DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      DiffusionMFD(plist, mesh),
      DiffusionWithGravity(mesh),
      Diffusion(mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
  }
  
  DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<Operator>& global_op,
                          const AmanziGeometry::Point& g) :
      DiffusionMFD(plist, global_op),
      DiffusionWithGravity(global_op),
      Diffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
  }
  
  DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const AmanziGeometry::Point& g) :
      DiffusionMFD(plist, mesh),
      DiffusionWithGravity(mesh),
      Diffusion(mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
  }
  
  DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<Operator>& global_op,
                          double rho, const AmanziGeometry::Point& g) :
      DiffusionMFD(plist, global_op),
      DiffusionWithGravity(global_op),
      Diffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
    SetDensity(rho);
  }

  DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          double rho, const AmanziGeometry::Point& g) :
      DiffusionMFD(plist, mesh),
      DiffusionWithGravity(mesh),
      Diffusion(mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
    SetDensity(rho);
  }

  // main members
  // -- required by the base class
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);
  // virtual void UpdateFluxNonManifold(const CompositeVector& u, CompositeVector& flux);

  // -- problem initialiation
  using DiffusionMFD::Setup;
  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             double rho, const AmanziGeometry::Point& g) {
    SetDensity(rho);
    SetGravity(g);
    DiffusionMFD::SetTensorCoefficient(K);
    DiffusionMFD::SetScalarCoefficient(k, dkdp);
  } 

  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const Teuchos::RCP<const CompositeVector>& rho,
             const AmanziGeometry::Point& g) {
    SetDensity(rho);
    SetGravity(g);
    DiffusionMFD::SetTensorCoefficient(K);
    DiffusionMFD::SetScalarCoefficient(k, dkdp);
  }

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const;

 protected:
  virtual void AddGravityToRHS_();
  inline AmanziGeometry::Point GravitySpecialDirection_(int f) const;
  void Init_(Teuchos::ParameterList& plist);

 protected:
  bool gravity_special_projection_;
  int gravity_method_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

