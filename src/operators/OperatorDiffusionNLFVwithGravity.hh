/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_NLFV_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_NVFV_WITH_GRAVITY_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Operators
#include "OperatorDiffusionNLFV.hh"
#include "OperatorDiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionNLFVwithGravity : public OperatorDiffusionNLFV,
                                         public OperatorDiffusionWithGravity {
 public:
  OperatorDiffusionNLFVwithGravity(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<Operator>& global_op,
                                   double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionNLFV(plist, global_op),
      OperatorDiffusionWithGravity(global_op),
      OperatorDiffusion(global_op)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV_GRAVITY;
    SetGravity(g);
    SetDensity(rho);
  }

  OperatorDiffusionNLFVwithGravity(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                   double rho, const AmanziGeometry::Point& g) :
      OperatorDiffusionNLFV(plist, mesh),
      OperatorDiffusionWithGravity(mesh),
      OperatorDiffusion(mesh)
  {
    operator_type_ = OPERATOR_DIFFUSION_NLFV_GRAVITY;
    SetGravity(g);
    SetDensity(rho);
  }

  // main virtual members 
  // -- setup
  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             double rho, const AmanziGeometry::Point& g) {
    SetGravity(g);
    SetDensity(rho);
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
  } 

  void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
             const Teuchos::RCP<const CompositeVector>& k,
             const Teuchos::RCP<const CompositeVector>& dkdp,
             const Teuchos::RCP<const CompositeVector>& rho,
             const AmanziGeometry::Point& g) {
    SetGravity(g);
    SetDensity(rho);
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
  }

  // -- create an operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  // -- modify an operator
  virtual void ModifyMatrices(const CompositeVector& u) {};
  virtual void ScaleMassMatrices(double s) {};

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const {};

  // virtual members from the base NLFV class
  // -- solution can be modified on boundary faces. This reflects specifics
  //    of nonlinear FV schemes.
  virtual double MapBoundaryValue_(int f, double u);
};

}  // namespace Operators
}  // namespace Amanzi

#endif
