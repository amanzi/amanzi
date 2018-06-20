/*
  This is the flow component of the Amanzi code.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "MatrixMFD_TPFA_ScaledConstraint.hh"

namespace Amanzi {
namespace Operators {


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD_TPFA_ScaledConstraint::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {

  if (Krel == Teuchos::null || !Krel->HasComponent("face")) {
    MatrixMFD_TPFA::CreateMFDstiffnessMatrices(Krel);
    Krel_->PutScalar(1.);

  } else {
    // tag global matrices as invalid
    MarkLocalMatricesAsChanged_();

    // store a copy of Krel on faces
    Krel->ScatterMasterToGhosted("face");
    *Krel_ = *(*Krel->ViewComponent("face", true))(0);

    int dim = mesh_->space_dimension();
    WhetStone::MFD3D_Diffusion mfd(mesh_);
    AmanziMesh::Entity_ID_List faces;

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    if (Aff_cells_.size() != ncells) {
      Aff_cells_.resize(static_cast<size_t>(ncells));
    }
    if (Afc_cells_.size() != ncells) {
      Afc_cells_.resize(static_cast<size_t>(ncells));
    }
    if (Acf_cells_.size() != ncells) {
      Acf_cells_.resize(static_cast<size_t>(ncells));
    }
    if (Acc_cells_.size() != ncells) {
      Acc_cells_.resize(static_cast<size_t>(ncells));
      Acc_ = Teuchos::rcp(new Epetra_Vector(View,mesh_->cell_map(false),&Acc_cells_[0]));
    }     

    for (int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseMatrix& Mff = Mff_cells_[c];
      Teuchos::SerialDenseMatrix<int, double> Bff(nfaces,nfaces);
      Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

      if (Krel->HasComponent("cell")) {
        const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);
        for (int n=0; n!=nfaces; ++n) {
          Bff(n, n) = Mff(n,n) * Krel_c[0][c];
        }
      } else {
        for (int n=0; n!=nfaces; ++n) {
          Bff(n, n) = Mff(n,n);
        }
      }

      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);
      double matsum = 0.0;  // elimination of mass matrix
      for (int n=0; n!=nfaces; ++n) {
        double rowsum = Bff(n, n);
        double colsum = Bff(n, n) * Krel_f[0][faces[n]];
        Bcf(n) = -colsum;
        Bfc(n) = -rowsum;
        matsum += colsum;
      }

      Aff_cells_[c] = Bff;
      Afc_cells_[c] = Bfc;
      Acf_cells_[c] = Bcf;
      Acc_cells_[c] = matsum;
    }
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi
