/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This test adds mode components to the nonlinear cofficient and overrides
  calculation of a local stiffness matrix.
*/

#ifndef AMANZI_OPERATOR_MY_PDE_DIFFUSION_MFD_HH_
#define AMANZI_OPERATOR_MY_PDE_DIFFUSION_MFD_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"

#include "PDE_DiffusionMFD.hh"

namespace Amanzi {
namespace Operators {

class MyPDE_DiffusionMFD : public virtual PDE_DiffusionMFD {
 public:
  MyPDE_DiffusionMFD(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(mesh), PDE_DiffusionMFD(plist, mesh){};

  // -- To calculate elemetal matrices, we can use input parameters flux
  //    and u from the previous nonlinear iteration. Otherwise, use null-pointers.
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;
};


/* ******************************************************************
* Second-order reconstruction of little k inside mesh cells.
* This member of DIVK-pamily of methods requires to recalcualte all
* mass matrices.
****************************************************************** */
void
MyPDE_DiffusionMFD::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                   const Teuchos::Ptr<const CompositeVector>& u)
{
  AMANZI_ASSERT(!scaled_constraint_);

  // preparing little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_grad = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_twin = Teuchos::null;
  if (k_ != Teuchos::null) {
    k_cell = k_->ViewComponent("cell");
    k_face = k_->ViewComponent("face", true);
    k_grad = k_->ViewComponent("grad");
    if (k_->HasComponent("twin")) k_twin = k_->ViewComponent("twin", true);
  }

  // update matrix blocks
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix Wff;

  AmanziMesh::Entity_ID_List cells;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c < ncells_owned; c++) {
    // mean value and gradient of nonlinear factor
    double kc = (*k_cell)[0][c];
    AmanziGeometry::Point kgrad(dim);
    for (int i = 0; i < dim; i++) kgrad[i] = (*k_grad)[i][c];

    // upwinded values of nonlinear factor
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();
    std::vector<double> kf(nfaces, 1.0);
    if (k_twin == Teuchos::null) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else {
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        kf[n] = (c == cells[0]) ? (*k_face)[0][f] : (*k_twin)[0][f];
      }
    }

    if (K_.get()) Kc = (*K_)[c];
    mfd.MassMatrixInverseDivKScaled(c, Kc, kc, kgrad, Wff);

    WhetStone::DenseMatrix Acell(nfaces + 1, nfaces + 1);

    double matsum = 0.0;
    for (int n = 0; n < nfaces; n++) {
      double rowsum = 0.0;
      for (int m = 0; m < nfaces; m++) {
        double tmp = Wff(n, m) * kf[n] * kf[m];
        rowsum += tmp;
        Acell(n, m) = tmp;
      }

      Acell(n, nfaces) = -rowsum;
      Acell(nfaces, n) = -rowsum;
      matsum += rowsum;
    }
    Acell(nfaces, nfaces) = matsum;
    local_op_->matrices[c] = Acell;
  }
}

} // namespace Operators
} // namespace Amanzi


#endif
