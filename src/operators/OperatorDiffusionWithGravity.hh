/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete gravity operator blended with the diffusion operator.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_WITH_GRAVITY_HH_

#include "Epetra_IntVector.h"

#include "tensor.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionMFD.hh"


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionWithGravity : public OperatorDiffusionMFD {
 public:
  OperatorDiffusionWithGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<Operator>& global_op) :
      OperatorDiffusionMFD(plist, global_op)
  {
    Init_(plist);
  }

  OperatorDiffusionWithGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      OperatorDiffusionMFD(plist, mesh)
  {
    Init_(plist);
  }

  OperatorDiffusionWithGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      OperatorDiffusionMFD(plist, mesh)
  {
    Init_(plist);
  }
  
  // main members
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  virtual void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }
  virtual void SetDensity(double rho) {
    scalar_rho_ = true;
    rho_ = rho;
  }
  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho) {
    scalar_rho_ = false;
    rho_cv_ = rho;
  }
  
 protected:
  virtual void AddGravityToRHS_();
  inline AmanziGeometry::Point GravitySpecialDirection_(int f) const;
  void Init_(Teuchos::ParameterList& plist);

 protected:
  AmanziGeometry::Point g_;
  bool gravity_special_projection_;
  int gravity_method_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

