/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_WITH_BND_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_NLFV_WITH_BND_WITH_GRAVITY_HH_

#include <strings.h>

// TPLs
#include "Ifpack.h" 
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseMatrix.hh"
#include "Preconditioner.hh"

// Operators
#include "PDE_DiffusionNLFVwithBndFaces.hh"
#include "PDE_DiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class PDE_DiffusionNLFVwithBndFacesGravity : public PDE_DiffusionNLFVwithBndFaces,
                                             public PDE_DiffusionWithGravity {
 public:
  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                                       const Teuchos::RCP<Operator>& global_op,
                                       double rho, const AmanziGeometry::Point& g) :
      PDE_Diffusion(global_op),
      PDE_DiffusionNLFVwithBndFaces(plist, global_op),
      PDE_DiffusionWithGravity(global_op)
  {
    pde_type_ = PDE_DIFFUSION_NLFVFACES_GRAVITY;
    SetGravity(g);
    SetDensity(rho);
  }

  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                                       const Teuchos::RCP<Operator>& global_op) :
      PDE_Diffusion(global_op),
      PDE_DiffusionNLFVwithBndFaces(plist, global_op),
      PDE_DiffusionWithGravity(global_op)
  {
    pde_type_ = PDE_DIFFUSION_NLFVFACES_GRAVITY;
  }  

  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                               double rho, const AmanziGeometry::Point& g) :
      PDE_Diffusion(mesh),
      PDE_DiffusionNLFVwithBndFaces(plist, mesh),
      PDE_DiffusionWithGravity(mesh)
  {
    pde_type_ = PDE_DIFFUSION_NLFVFACES_GRAVITY;
    SetGravity(g);
    SetDensity(rho);
  }

  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Diffusion(mesh),
      PDE_DiffusionNLFVwithBndFaces(plist, mesh),
      PDE_DiffusionWithGravity(mesh)
  {
    pde_type_ = PDE_DIFFUSION_NLFVFACES_GRAVITY;
  }

  PDE_DiffusionNLFVwithBndFacesGravity(Teuchos::ParameterList& plist,
                                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                       const Teuchos::RCP<const CompositeVector>& rho,
                                       const AmanziGeometry::Point& g) :
      PDE_Diffusion(mesh),
      PDE_DiffusionNLFVwithBndFaces(plist, mesh),
      PDE_DiffusionWithGravity(mesh)
  {
    pde_type_ = PDE_DIFFUSION_NLFVFACES_GRAVITY;
    SetGravity(g);
    SetDensity(rho);
  }

  // main virtual members 
  // -- setup
  using PDE_Diffusion::Setup;
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
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  // -- after solving the problem: postrocessing
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;
  virtual void UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
                                     const Teuchos::Ptr<CompositeVector>& flux) override {};

  // -- modify an operator
  virtual void ModifyMatrices(const CompositeVector& u) override {};
  virtual void ScaleMassMatrices(double s) override {};

  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  virtual double ComputeGravityFlux(int f) const override {
    AMANZI_ASSERT(0);
    return 0.;
  };

  // virtual members from the base NLFV class
  // -- solution can be modified on boundary faces. This reflects specifics
  //    of nonlinear FV schemes.
  // virtual double MapBoundaryValue_(int f, double u) override;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
