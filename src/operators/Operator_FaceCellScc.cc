/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are CELL + FACE, but which assembles the
  CELL only system and Schur complements the face.

  This uses special assembly. Apply is done as if we had the full 
  FACE+CELL system. SymbolicAssembly() is done as if we had the CELL 
  system, but with an additional step to get the layout due to the 
  Schur'd system on FACE+CELL. Assemble, however, is done using a 
  totally different approach.
*/

#include <vector>

#include "EpetraExt_RowMatrixOut.h"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Preconditioner.hh"

#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "Op.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"
#include "Operator_FaceCellScc.hh"

namespace Amanzi {
namespace Operators {

// Special Apply Inverse required to deal with schur complement
int Operator_FaceCellScc::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  int ierr(0);
  Y.PutScalarGhosted(0.0);

  // apply preconditioner inversion
  // const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  // Temporary cell and face vectors.
  CompositeVector T(X);
  Epetra_MultiVector& Tf = *T.ViewComponent("face", true);
  Epetra_MultiVector& Tc = *T.ViewComponent("cell");

  // lump Aff
  AmanziMesh::Entity_ID_List faces;

  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL |
                               OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      for (int c = 0; c != ncells_owned; c++) {
        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        WhetStone::DenseMatrix& Acell = (*it)->matrices[c];

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          Tf[0][f] += Acell(n, n);
        }
      }
    }
  }
  T.GatherGhostedToMaster("face");

  // FORWARD ELIMINATION:  Tc = Xc - Acf inv(Aff) Xf
  T.ScatterMasterToGhosted("face");
  X.ScatterMasterToGhosted("face");

  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL |
                               OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      for (int c = 0; c < ncells_owned; c++) {
        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        WhetStone::DenseMatrix& Acell = (*it)->matrices[c];

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          Tc[0][c] -= Acell(n, nfaces) / Tf[0][f] * Xf[0][f];
        }
      }
    }
  }

  // Solve the Schur complement system Yc = inv(Scc) * Tc.
  {
    preconditioner_->ApplyInverse(T, Y);
  }

  // BACKWARD SUBSTITUTION:  Yf = inv(Aff) (Xf - Afc Yc)
  {
    const Epetra_MultiVector& Yc = *Y.ViewComponent("cell");
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    Yf = Xf;

    for (const_op_iterator it = begin(); it != end(); ++it) {
      if ((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL |
                                 OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
        for (int c = 0; c < ncells_owned; c++) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          WhetStone::DenseMatrix& Acell = (*it)->matrices[c];

          for (int n = 0; n < nfaces; n++) {
            int f = faces[n];
            Yf[0][f] -= Acell(nfaces, n) * Yc[0][c];
          }
        }
      }
    }
  }
  Y.GatherGhostedToMaster("face", Add);

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    for (int f = 0; f < nfaces_owned; f++) Yf[0][f] /= Tf[0][f];
  }

  return ierr;
}


// Special AssembleMatrix required to deal with schur complement
void Operator_FaceCellScc::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
                                          int my_block_row, int my_block_col) const
{
  // first check preconditions -- Scc must have exactly one face-based schema (a FACE+CELL)
  int num_with_faces = 0;
  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ & OPERATOR_SCHEMA_DOFS_FACE) {
      num_with_faces++;
    }
  }
  if (num_with_faces == 0 || num_with_faces > 1) {
    Errors::Message msg("Schur complement to Scc must have exactly one face-based schema for diffusion.");
    Exceptions::amanzi_throw(msg);
  }

  
  // schur complement can be split into a CELL_CELL schema, Acc, and a
  // FACE_CELL schema, -Acf*(Aff^-1)*Afc.
  int i_schur = 0;
  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL |
                               OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      AMANZI_ASSERT((*it)->matrices.size() == ncells_owned);

      // create or get extra ops, and keep them for future use
      Teuchos::RCP<Op_Cell_Cell> diag_op;
      Teuchos::RCP<Op_Face_Cell> schur_op;

      if (diag_ops_.size() > i_schur) {
        // get existing
        diag_op = diag_ops_[i_schur];
        schur_ops_[i_schur];
      } else {
        // create and fill
        AMANZI_ASSERT(i_schur == diag_ops_.size());
        std::string name = "Scc alt as CELL_CELL";
        diag_op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_));
        diag_ops_.push_back(diag_op);

        name = "Scc alt as FACE_CELL";
        schur_op = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
        schur_ops_.push_back(schur_op);
        AmanziMesh::Entity_ID_List cells;
        for (int f = 0; f != nfaces_owned; ++f) {
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
          int ncells = cells.size();
          schur_op->matrices[f] = WhetStone::DenseMatrix(ncells, ncells);
        }
      }

      // set the size of the local val vectors
      Epetra_MultiVector& diag = *diag_op->diag;

      // populate the diagonal component
      for (int c = 0; c != ncells_owned; ++c) {
        WhetStone::DenseMatrix& Acell = (*it)->matrices[c];
        int n = Acell.NumCols() - 1;
        diag[0][c] = (*it)->matrices[c](n,n);
        if (diag[0][c] > 1.e30 || diag[0][c] < -1.e30) {
          AMANZI_ASSERT(0);
        }
      }

      // populate the schur component
      // -- populate coefficients that form the transmissibilities
      const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

      CompositeVectorSpace cv_space;
      cv_space.SetMesh(mesh_);
      cv_space.SetGhosted(true);
      cv_space.SetComponent("face", AmanziMesh::FACE, 2);

      CompositeVector T(cv_space, true);
      Epetra_MultiVector& Ttmp = *T.ViewComponent("face", true);

      AmanziMesh::Entity_ID_List cells, faces;

      Ttmp.PutScalar(0.0);
      for (int c = 0; c < ncells_owned; c++) {
        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        int c0 = cmap_wghost.GID(c);
        WhetStone::DenseMatrix& Acell = (*it)->matrices[c];
        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

          if (cells.size() == 1) {
            Ttmp[0][f] = Acell(n, nfaces);
          } else {
            int c1 = cmap_wghost.GID(cells[0]);
            int c2 = cmap_wghost.GID(cells[1]);
            int i = (c0 == std::min(c1, c2)) ? 0 : 1;
            Ttmp[i][f] = Acell(n, nfaces);
          }
        }
      }
      T.GatherGhostedToMaster();

      // -- populate the local matrices
      int gid[2];
      double a1, a2;

      std::vector<WhetStone::DenseMatrix>& mats = schur_op->matrices;

      for (int f = 0; f < nfaces_owned; f++) {
        WhetStone::DenseMatrix& mat = mats[f];

        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n=0; n!=ncells; ++n) {
          gid[n] = cmap_wghost.GID(cells[n]);
        }

        a1 = Ttmp[0][f];
        a2 = Ttmp[1][f];
        if (ncells == 2 && gid[0] > gid[1]) {
          a1 = Ttmp[1][f];
          a2 = Ttmp[0][f];
        }

        double coef = fabs(a1 + a2);
        if (coef == 0.0) continue;

        mat(0,0) = -a1 * a1 / coef;
        if (mat(0,0) > 1.e30 || mat(0,0) < -1.e30) {
          AMANZI_ASSERT(0);
        }
        if (ncells > 1) {
          mat(0,1) = -a1 * a2 / coef;
          mat(1,0) = -a1 * a2 / coef;
          mat(1,1) = -a2 * a2 / coef;
          if (mat(1,0) > 1.e30 || mat(1,0) < -1.e30
              || mat(0,1) > 1.e30 || mat(0,1) < -1.e30
              || mat(1,1) > 1.e30 || mat(1,1) < -1.e30) {
            AMANZI_ASSERT(0);
          }
        }
      }

      // assemble
      diag_op->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
      AMANZI_ASSERT(schur_op->matrices.size() == nfaces_owned);
      schur_op->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
      i_schur++;
    } else {
      (*it)->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
    }
  } // for
}


// visit method for Apply -- this is identical to Operator_FaceCell's
// version.
int Operator_FaceCellScc::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                            const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  Y.PutScalarGhosted(0.);
  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

    AmanziMesh::Entity_ID_List faces;
    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
      for (int n=0; n!=nfaces; ++n) {
        v(n) = Xf[0][faces[n]];
      }
      v(nfaces) = Xc[0][c];

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n=0; n!=nfaces; ++n) {
        Yf[0][faces[n]] += av(n);
      }
      Yc[0][c] += av(nfaces);
    }
  }
  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
 * Create a global matrix.
 ****************************************************************** */
void Operator_FaceCellScc::SymbolicAssembleMatrix()
{
  // SuperMap for Sff is face only
  CompositeVectorSpace smap_space;
  smap_space.SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  smap_ = createSuperMap(smap_space);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema(), 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
          smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  Operator::SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  AMANZI_ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


// visit method for sparsity structure of Schur complement
void Operator_FaceCellScc::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                                    const SuperMap& map, GraphFE& graph,
                                                    int my_block_row, int my_block_col) const
{
  std::string name = "Scc alt as FACE_CELL";
  Op_Face_Cell schur_op(name, mesh_);
  Operator_Cell::SymbolicAssembleMatrixOp(schur_op, map, graph, my_block_row, my_block_col);

  name = "Scc alt as CELL_CELL";
  Op_Cell_Cell diag_op(name, mesh_);
  Operator_Cell::SymbolicAssembleMatrixOp(diag_op, map, graph, my_block_row, my_block_col);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator> Operator_FaceCellScc::Clone() const {
  return Teuchos::rcp(new Operator_FaceCellScc(*this));
}

}  // namespace Operators
}  // namespace Amanzi

