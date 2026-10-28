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

#include "Kokkos_Core.hpp"

#include "DenseMatrix.hh"
#include "MeshCache.hh"
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
GravitySpecialDirection(const AmanziMesh::MeshCache& mesh, int f)
{
  auto cells = mesh.getFaceCells(f);
  int ncells = cells.size();

  if (ncells == 2) {
    return mesh.getCellCentroid(cells[1]) - mesh.getCellCentroid(cells[0]);
  } else {
    return mesh.getFaceCentroid(f) - mesh.getCellCentroid(cells[0]);
  }
}
} // namespace Impl

class BCs;

class PDE_DiffusionMFDwithGravity : public PDE_DiffusionMFD, public PDE_DiffusionWithGravity {
 public:
  PDE_DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                              const Teuchos::RCP<Operator>& global_op)
    : PDE_Diffusion(plist, global_op),
      PDE_DiffusionMFD(plist, global_op),
      PDE_DiffusionWithGravity(plist, global_op)
  {}

  PDE_DiffusionMFDwithGravity(Teuchos::ParameterList& plist,
                              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(plist, mesh),
      PDE_DiffusionMFD(plist, mesh),
      PDE_DiffusionWithGravity(plist, mesh)
  {}

  // delete this to avoid accidental capture in lambdas
  PDE_DiffusionMFDwithGravity(const PDE_DiffusionMFDwithGravity& other) = delete;

  virtual void Init() override;

  // main members
  // -- required by the base class
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;


  virtual void AddGravityToRHS_();

  // Lightweight, KOKKOS_FUNCTION-callable, on-demand per-face
  // {transmissibility, gravity flux} lookup. Transmissibility looks up the
  // DIAGONAL entry of the already-cached per-cell dense mass matrix
  // Wff_cells_ -- exact in TPFA submode (where Wff is diagonal by
  // construction), an acknowledged TPFA-consistent APPROXIMATION otherwise
  // (general MFD's flux through a face genuinely depends on all faces of
  // the cell via the full dense Wff, not just the diagonal); callers using
  // this on a general-MFD operator must treat it as such at the specific
  // face(s) where they use it, not as a replacement for MFD's own interior
  // flux computation (untouched). Gravity flux is a genuine closed-form
  // per-(cell,face) computation (mirrors AddGravityToRHS_/UpdateFlux's
  // fv_flag branch, specialized to the single-neighbor-cell boundary-face
  // case), WITHOUT the mobility (kr_n) factor -- matching the convention of
  // FV's gravity_term_, the caller applies mobility itself. Meant to be
  // captured by value into a device loop (see
  // WithTransmissibilityGravityFunctor in PDE_DiffusionFactory.hh), not
  // called through a virtual pointer per face.
  struct TransGravFunctor {
    DenseMatrix_Vector Wff_cells;
    AmanziMesh::MeshCache mesh;
    bool have_K;
    TensorVector K;
    CompositeVector::cView_type rho_c;
    AmanziGeometry::Point g;

    KOKKOS_INLINE_FUNCTION Kokkos::pair<double, double> operator()(int f) const
    {
      // A coupling face is a boundary face -- exactly one neighbor cell.
      auto cells = mesh.getFaceCells(f);
      int c = cells(0);
      auto cfaces = mesh.getCellFaces(c);
      int n = -1;
      for (int i = 0; i != cfaces.extent(0); ++i) {
        if (cfaces(i) == f) {
          n = i;
          break;
        }
      }
      auto Wff = Wff_cells[c];
      double trans = Wff(n, n);

      const AmanziGeometry::Point& normal = mesh.getFaceNormal(f);
      AmanziGeometry::Point Kcg = have_K ? K.at(c) * g : g;
      double grav = (Kcg * normal) * rho_c(c, 0);
      return { trans, grav };
    }
  };

  TransGravFunctor GetTransmissibilityGravityFunctor() const
  {
    if (!mass_matrices_initialized_) {
      Errors::Message msg(
        "PDE_DiffusionMFDwithGravity::GetTransmissibilityGravityFunctor: mass matrices have not "
        "yet been computed -- call UpdateMatrices() first.");
      Exceptions::amanzi_throw(msg);
    }
    static const TensorVector kEmptyK;
    bool have_K = K_ != Teuchos::null;
    return { Wff_cells_, mesh_->getCache(), have_K, have_K ? *K_ : kEmptyK, DensityCells(true), g_ };
  }

 protected:
  bool gravity_special_projection_;
  int gravity_method_;
};

} // namespace Operators
} // namespace Amanzi

#endif
