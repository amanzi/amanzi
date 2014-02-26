/*
  This is the Operator component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"

#include "DenseVector.hh"
#include "PreconditionerFactory.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Default constructor.
****************************************************************** */
Operator::Operator(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) 
    : cvs_(cvs), data_validity_(true) 
{
  mesh_ = cvs_->Mesh();
  rhs_ = Teuchos::rcp(new CompositeVector(*cvs_, true));
  diagonal_block_ = Teuchos::rcp(new CompositeVector(*cvs_, true));
  diagonal_block_->PutScalar(0.0);

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
};


/* ******************************************************************
* Second constructor.
****************************************************************** */
Operator::Operator(const Operator& op) 
    : data_validity_(true), 
      mesh_(op.mesh_),
      cvs_(op.cvs_), 
      matrix_blocks_(op.matrix_blocks_), 
      matrix_blocks_type_(op.matrix_blocks_type_), 
      diagonal_block_(op.diagonal_block_),
      rhs_(op.rhs_),
      A_(op.A_),
      preconditioner_(op.preconditioner_)
{
  op.data_validity_ = false;
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
}


/* ******************************************************************
* Initialization of the operator.                                           
****************************************************************** */
void Operator::Init()
{
  diagonal_block_->PutScalar(0.0);
  rhs_->PutScalar(0.0);
}


/* ******************************************************************
* Create a face-based global matrix.
****************************************************************** */
void Operator::SymbolicAssembleFaces()
{
  // map used with the matrix
  const Epetra_Map& map = mesh_->face_map(false);
  const Epetra_Map& map_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
  Epetra_FECrsGraph ff_graph(Copy, map, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int faces_LID[OPERATOR_MAX_FACES];
  int faces_GID[OPERATOR_MAX_FACES];

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      faces_LID[n] = faces[n];
      faces_GID[n] = map_wghost.GID(faces_LID[n]);
    }
    ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }
  ff_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  A_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  A_->GlobalAssemble();
}


/* ******************************************************************
* Assemble elemental face-based matrices into four global matrices. 
****************************************************************** */
void Operator::AssembleStencilMFD_Faces()
{
  // find location of face-based matrices
  int m, nblocks = matrix_blocks_.size();
  for (int n = 0; n < nblocks; n++) {
    if (matrix_blocks_type_[n] == OPERATOR_STENCIL_CELL_CF_CF) {
      m = n;
      break;
    }
  }

  std::vector<WhetStone::DenseMatrix>& matrix = *matrix_blocks_[m];
  Epetra_MultiVector& diagonal = *diagonal_block_->ViewComponent("face");

  // populate the matrix
  A_->PutScalar(0.0);

  const Epetra_Map& map = mesh_->face_map(false);
  const Epetra_Map& map_wghost = mesh_->face_map(true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int faces_LID[OPERATOR_MAX_FACES];
  int faces_GID[OPERATOR_MAX_FACES];

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      faces_LID[n] = faces[n];
      faces_GID[n] = map_wghost.GID(faces_LID[n]);
    }
    A_->SumIntoGlobalValues(nfaces, faces_GID, matrix[c].Values());
  }
  A_->GlobalAssemble();

  // Add diagonal
  diagonal_block_->GatherGhostedToMaster("face", Add);

  Epetra_Vector tmp(map);
  A_->ExtractDiagonalCopy(tmp);
  tmp.Update(1.0, diagonal, 1.0);
  A_->ReplaceDiagonalValues(tmp);

  // Assemble all right-hand sides
  rhs_->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Applies boundary conditions to matrix_blocks and update the
* right-hand side and the diagonal block.                                           
****************************************************************** */
void Operator::ApplyBCs(std::vector<int>& bc_model, std::vector<double>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int nblocks = matrix_blocks_.size();
  for (int n = 0; n < nblocks; n++) {
    std::vector<WhetStone::DenseMatrix>& matrix = *matrix_blocks_[n];
    int type = matrix_blocks_type_[n];

    if (type == OPERATOR_STENCIL_CELL_CF_CF) {
      Epetra_MultiVector& rhs_face = *rhs_->ViewComponent("face");
      Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");
      Epetra_MultiVector& diagonal = *diagonal_block_->ViewComponent("face");

      for (int c = 0; c < ncells_owned; c++) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        WhetStone::DenseMatrix& Acell = matrix[c];

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          double value = bc_values[f];

          if (bc_model[f] == OPERATOR_BC_FACE_DIRICHLET) {
            for (int m = 0; m < nfaces; m++) {
              rhs_face[0][faces[m]] -= Acell(m, n) * value;
              Acell(n, m) = Acell(m, n) = 0.0;
            }
            rhs_face[0][f] = value;
            diagonal[0][f] = 1.0;

            rhs_cell[0][c] -= Acell(nfaces, n) * value;
            Acell(nfaces, n) = 0.0;
            Acell(n, nfaces) = 0.0;
          }
        }
      }
    }
  }
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
******************************************************************* */
int Operator::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  // Multiply by the diagonal block.
  Y.Multiply(1.0, *diagonal_block_, X, 0.0);

  // Multiply by the remaining matrix blocks.
  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

  int nblocks = matrix_blocks_.size();
  for (int n = 0; n < nblocks; n++) {
    std::vector<WhetStone::DenseMatrix>& matrix = *matrix_blocks_[n];
    int type = matrix_blocks_type_[n];

    if (type == OPERATOR_STENCIL_CELL_CF_CF) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

      Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
      Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

      for (int c = 0; c < ncells_owned; c++) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
        for (int n = 0; n < nfaces; n++) {
          v(n) = Xf[0][faces[n]];
        }
        v(nfaces) = Xc[0][c];

        WhetStone::DenseMatrix& Acell = matrix[c];
        Acell.Multiply(v, av, false);

        for (int n = 0; n < nfaces; n++) {
          Yf[0][faces[n]] += av(n);
        }
        Yc[0][c] += av(nfaces);
      } 
    } else if (type == OPERATOR_STENCIL_FACE_C_C) {
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);
      Epetra_MultiVector& Yc = *Y.ViewComponent("cell", true);

      for (int f = 0; f < nfaces_owned; f++) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        int ncells = cells.size();

        WhetStone::DenseVector v(ncells), av(ncells);
        for (int n = 0; n < ncells; n++) {
          v(n) = Xc[0][cells[n]];
        }

        WhetStone::DenseMatrix& Aface = matrix[f];
        Aface.Multiply(v, av, false);

        for (int n = 0; n < ncells; n++) {
          Yc[0][cells[n]] += av(n);
        }
      } 
    }
  }
  Y.GatherGhostedToMaster(Add);

  return 0;
}


/* ******************************************************************
 * * Initialization of the preconditioner                                                 
 * ****************************************************************** */
void Operator::InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, plist);
  preconditioner_->Update(A_);
}

}  // namespace Operators
}  // namespace Amanzi



