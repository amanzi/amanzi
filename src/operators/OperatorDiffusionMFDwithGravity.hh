/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete gravity operator blended with the MFD diffusion operator.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_MFD_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_MFD_WITH_GRAVITY_HH_

#include "Epetra_IntVector.h"

#include "tensor.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionMFD.hh"


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionMFDwithGravity : public OperatorDiffusionMFD {
 public:
  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<Operator>& global_op,
                                  double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, global_op),
      rho_(rho),
      g_(g)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
  }

  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                  double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, mesh),
      rho_(rho),
      g_(g)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
  }

  OperatorDiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                  Teuchos::RCP<const CompositeVector>& rho_cv,
                                  const AmanziGeometry::Point& g) :
      OperatorDiffusionMFD(plist, mesh),
      rho_cv_(rho_cv),
      g_(g)
  {
    operator_type_ = OPERATOR_DIFFUSION_MFD_GRAVITY;
    Init_(plist);
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
    rho_ = rho;
    g_  = g;
    OperatorDiffusionMFD::Setup(K);
    Setup(k, dkdp);
  } 

  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const Teuchos::RCP<const CompositeVector>& rho,
             const AmanziGeometry::Point& g) {
    rho_cv_ = rho;
    g_  = g;
    Setup(K);
    Setup(k, dkdp);
  }

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const;

 protected:
  virtual void AddGravityToRHS_();
  inline AmanziGeometry::Point GravitySpecialDirection_(int f) const;
  void Init_(Teuchos::ParameterList& plist);

 protected:
  double rho_;
  Teuchos::RCP<const CompositeVector> rho_cv_;

  AmanziGeometry::Point g_;
  bool gravity_special_projection_;
  int gravity_method_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

