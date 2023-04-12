// PDE_DiffusionMFDwithGravity: Discrete gravity operator blended with the MFD diffusion operator.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_MFD_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_MFD_WITH_GRAVITY_HH_

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionWithGravity.hh"

/*!
Additional options for MFD with the gravity term include:

* `"gravity term discretization`" ``[string]`` selects a model for discretizing the
  gravity term. Available options are `"hydraulic head`" (default) and `"finite volume`".
  The first option starts with equation for the shifted solution, i.e. the hydraulic
  head, and derives gravity discretization by the reserve shifting.
  The second option is based on the divergence formula.
*/

namespace Amanzi {
namespace Operators {

namespace Impl {
  // Developments
  // -- interface to solvers for treating nonlinear BCs.
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  GravitySpecialDirection(const AmanziMesh::Mesh* const mesh, int f)
  {
    auto cells = mesh->getFaceCells(f);
    int ncells = cells.size();

    if (ncells == 2) {
      return mesh->getCellCentroid(cells[1]) - mesh->getCellCentroid(cells[0]);
    } else {
      return mesh->getFaceCentroid(f) - mesh->getCellCentroid(cells[0]);
    }
  }
} // namespace Impl

class BCs;

class PDE_DiffusionMFDwithGravity : public PDE_DiffusionMFD, public PDE_DiffusionWithGravity {

 public:
  PDE_DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                              const Teuchos::RCP<Operator>& global_op) :
    PDE_Diffusion(plist, global_op),
    PDE_DiffusionMFD(plist, global_op),
    PDE_DiffusionWithGravity(plist, global_op)
  {}

  PDE_DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    PDE_Diffusion(plist, mesh),
    PDE_DiffusionMFD(plist, mesh),
    PDE_DiffusionWithGravity(plist, mesh)
  {}

  virtual void Init() override;

  // main members
  // -- required by the base class
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

  // virtual void UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
  //                                    const Teuchos::Ptr<CompositeVector>& flux) override;


  virtual void AddGravityToRHS_();

 protected:
  bool gravity_special_projection_;
  int gravity_method_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

