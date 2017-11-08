// OperatorDiffusionMFDwithGravity: Discrete gravity operator blended with the MFD diffusion operator.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_MFD_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_MFD_WITH_GRAVITY_HH_

#include "Epetra_IntVector.h"

#include "Tensor.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionMFD.hh"
#include "OperatorDiffusionWithGravity.hh"

/*!
Additional options for MFD with the gravity term include:
  
* `"gravity term discretization`" [string] selects a model for discretizing the 
   gravity term. Available options are `"hydraulic head`" [default] and `"finite volume`". 
   The first option starts with equation for the shifted solution, i.e. the hydraulic head,
   and derives gravity discretization by the reserve shifting.
   The second option is based on the divergence formula.
*/


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionMFDwithGravity : public OperatorDiffusionMFD,
                                        public OperatorDiffusionWithGravity {
 public:
  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusionMFD(plist, global_op),
      OperatorDiffusionWithGravity(global_op),
      OperatorDiffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
  }
  
  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusionMFD(plist, mesh),
      OperatorDiffusionWithGravity(mesh),
      OperatorDiffusion(mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
  }
  
  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<Operator>& global_op,
                                  const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, global_op),
      OperatorDiffusionWithGravity(global_op),
      OperatorDiffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
  }
  
  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                  const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, mesh),
      OperatorDiffusionWithGravity(mesh),
      OperatorDiffusion(mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
  }
  
  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<Operator>& global_op,
                                  double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, global_op),
      OperatorDiffusionWithGravity(global_op),
      OperatorDiffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
    SetGravity(g);
    SetDensity(rho);
  }

  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                  double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, mesh),
      OperatorDiffusionWithGravity(mesh),
      OperatorDiffusion(mesh)
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

  // -- problem initialiation
  using OperatorDiffusionMFD::Setup;
  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             double rho, const AmanziGeometry::Point& g) {
    SetDensity(rho);
    SetGravity(g);
    OperatorDiffusionMFD::SetTensorCoefficient(K);
    OperatorDiffusionMFD::SetScalarCoefficient(k, dkdp);
  } 

  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const Teuchos::RCP<const CompositeVector>& rho,
             const AmanziGeometry::Point& g) {
    SetDensity(rho);
    SetGravity(g);
    OperatorDiffusionMFD::SetTensorCoefficient(K);
    OperatorDiffusionMFD::SetScalarCoefficient(k, dkdp);
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

