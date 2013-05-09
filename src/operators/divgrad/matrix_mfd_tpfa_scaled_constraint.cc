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

#include "matrix_mfd_tpfa_scaled_constraint.hh"

namespace Amanzi {
namespace Operators {


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD_TPFA_ScaledConstraint::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {

  if (Krel == Teuchos::null || !Krel->has_component("face")) {
    MatrixMFD_TPFA::CreateMFDstiffnessMatrices(Krel);
  } else {
    // tag global matrices as invalid
    assembled_schur_ = false;
    assembled_operator_ = false;

    // store a copy of Krel on faces
    *Krel_ = *(*Krel->ViewComponent("face", true))(0);

    int dim = mesh_->space_dimension();
    WhetStone::MFD3D mfd(mesh_);
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    Aff_cells_.clear();
    Afc_cells_.clear();
    Acf_cells_.clear();
    Acc_cells_.clear();

    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      Teuchos::SerialDenseMatrix<int, double>& Mff = Mff_cells_[c];
      Teuchos::SerialDenseMatrix<int, double> Bff(nfaces,nfaces);
      Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

      if (Krel->has_component("cell")) {
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

      Aff_cells_.push_back(Bff);
      Afc_cells_.push_back(Bfc);
      Acf_cells_.push_back(Bcf);
      Acc_cells_.push_back(matsum);
    }
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi
