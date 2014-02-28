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
  diagonal_ = Teuchos::rcp(new CompositeVector(*cvs_, true));
  diagonal_->PutScalar(0.0);

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
}


/* ******************************************************************
* Second constructor.
****************************************************************** */
Operator::Operator(const Operator& op) 
    : data_validity_(true), 
      mesh_(op.mesh_),
      cvs_(op.cvs_), 
      blocks_(op.blocks_), 
      blocks_schema_(op.blocks_schema_), 
      diagonal_(op.diagonal_),
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
  diagonal_->PutScalar(0.0);
  rhs_->PutScalar(0.0);
}


/* ******************************************************************
* Create a global matrix.
****************************************************************** */
void Operator::SymbolicAssembleMatrix(int schema)
{
  // create global Epetra map
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& vmap = mesh_->node_map(false);

  int ndof(0), ndof_global(0), offset(0);

  if (schema & OPERATOR_DOFS_FACE) ndof += nfaces_owned;
  if (schema & OPERATOR_DOFS_CELL) ndof += ncells_owned;
  if (schema & OPERATOR_DOFS_NODE) ndof += nnodes_owned;

  int* gids = new int[ndof];

  offset_global_[0] = 0;
  offset_my_[0] = 0;
  if (schema & OPERATOR_DOFS_FACE) {
    fmap.MyGlobalElements(&(gids[0]));
    offset += nfaces_owned;
    ndof_global += fmap.NumGlobalElements();
  } 

  offset_global_[1] = ndof_global;
  offset_my_[1] = offset;
  if (schema & OPERATOR_DOFS_CELL) {
    cmap.MyGlobalElements(&(gids[offset]));
    for (int c = 0; c < ncells_owned; c++) gids[offset + c] += ndof_global;
    offset += ncells_owned;
    ndof_global += cmap.NumGlobalElements();
  }

  offset_global_[2] = ndof_global;
  offset_my_[2] = offset;
  if (schema & OPERATOR_DOFS_NODE) {
    vmap.MyGlobalElements(&(gids[offset]));
    for (int v = 0; v < ncells_owned; v++) gids[offset + v] += ndof_global;
    offset += nnodes_owned;
    ndof_global += vmap.NumGlobalElements();
  }

  Epetra_Map* map = new Epetra_Map(ndof_global, ndof, gids, 0, cmap.Comm());
  delete [] gids;

  // estimate size of the matrix graph
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& vmap_wghost = mesh_->node_map(true);

  int row_size(0), dim = mesh_->space_dimension();
  if (schema & OPERATOR_DOFS_FACE) {
    int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
    row_size += 2 * i;
  }
  if (schema & OPERATOR_DOFS_CELL) {
    int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
    row_size += i + 1;
  }
  if (schema & OPERATOR_DOFS_NODE) {
    int i = (dim == 2) ? OPERATOR_QUAD_NODES : OPERATOR_HEX_NODES;
    row_size += 8 * i;
  }
    
  Epetra_FECrsGraph ff_graph(Copy, *map, row_size);
  delete map;

  // populate matrix graph
  AmanziMesh::Entity_ID_List cells, faces, nodes;
  std::vector<int> dirs;

  int lid[OPERATOR_MAX_NODES];
  int gid[OPERATOR_MAX_NODES];

  if (schema & OPERATOR_DOFS_FACE) {
    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; n++) {
        lid[n] = faces[n];
        gid[n] = fmap_wghost.GID(lid[n]);
      }
      ff_graph.InsertGlobalIndices(nfaces, lid, nfaces, gid);
    }
  }
  if (schema & OPERATOR_DOFS_CELL) {
    for (int f = 0; f < nfaces_owned; f++) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      for (int n = 0; n < ncells; n++) {
        lid[n] = offset_my_[1] + cells[n];
        gid[n] = offset_global_[1] + cmap_wghost.GID(cells[n]);
      }
      ff_graph.InsertGlobalIndices(ncells, lid, ncells, gid);
    }
  }
  if (schema & OPERATOR_DOFS_NODE) {
    for (int v = 0; v < nnodes_owned; v++) {
      mesh_->cell_get_nodes(v, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n < nnodes; n++) {
        lid[n] = offset_my_[2] + nodes[n];
        gid[n] = offset_global_[2] + vmap_wghost.GID(nodes[n]);
      }
      ff_graph.InsertGlobalIndices(nnodes, lid, nnodes, gid);
    }
  }

  ff_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrix
  A_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  A_->GlobalAssemble();
}


/* ******************************************************************
* Assemble elemental face-based matrices into four global matrices. 
****************************************************************** */
void Operator::AssembleMatrix(int schema)
{
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& vmap_wghost = mesh_->node_map(true);

  // populate matrix
  A_->PutScalar(0.0);

  AmanziMesh::Entity_ID_List faces, nodes;
  std::vector<int> dirs;

  int lid[OPERATOR_MAX_NODES + 1];
  int gid[OPERATOR_MAX_NODES + 1];

  int nblocks = blocks_.size();
  for (int n = 0; n < nblocks; n++) {
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[n];
    if ((blocks_schema_[n] & schema) != blocks_schema_[n]) {
      std::cout << "Schema " << schema 
                << " differs from the block schema " << blocks_schema_[n] << std::endl;
      ASSERT(0);
    }

    if (blocks_schema_[n] & OPERATOR_DOFS_FACE && 
        blocks_schema_[n] & OPERATOR_DOFS_CELL) {
      for (int c = 0; c < ncells_owned; c++) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int n = 0; n < nfaces; n++) {
          lid[n] = faces[n];
          gid[n] = fmap_wghost.GID(lid[n]);
        }
        lid[nfaces] = offset_my_[1] + c;
        gid[nfaces] = offset_global_[1] + cmap_wghost.GID(c);

        A_->SumIntoGlobalValues(nfaces + 1, gid, matrix[c].Values());
      }
    }

    if (blocks_schema_[n] & OPERATOR_DOFS_NODE) {
      for (int c = 0; c < ncells_owned; c++) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int n = 0; n < nnodes; n++) {
          lid[n] = offset_my_[2] + nodes[n];
          gid[n] = offset_global_[2] + vmap_wghost.GID(nodes[n]);
        }

        A_->SumIntoGlobalValues(nnodes, gid, matrix[c].Values());
      }
    }
  }

  A_->GlobalAssemble();

  // Add diagonal
  diagonal_->GatherGhostedToMaster(Add);
  Epetra_MultiVector& diag = *diagonal_->ViewComponent("face");

  Epetra_Vector tmp(A_->RowMap());
  A_->ExtractDiagonalCopy(tmp);
  tmp.Update(1.0, diag, 1.0);
  A_->ReplaceDiagonalValues(tmp);

  // Assemble all right-hand sides
  rhs_->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Applies boundary conditions to matrix_blocks and update the
* right-hand side and the diagonal block.                                           
****************************************************************** */
void Operator::ApplyBCs(std::vector<int>& bc_model, std::vector<double>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int nblocks = blocks_.size();
  for (int n = 0; n < nblocks; n++) {
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[n];
    int schema = blocks_schema_[n];

    if (schema & OPERATOR_DOFS_FACE) {
      Epetra_MultiVector& rhs_face = *rhs_->ViewComponent("face");
      Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");
      Epetra_MultiVector& diag = *diagonal_->ViewComponent("face");

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
            diag[0][f] = 1.0;

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
  Y.Multiply(1.0, *diagonal_, X, 0.0);

  // Multiply by the remaining matrix blocks.
  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

  int nblocks = blocks_.size();
  for (int n = 0; n < nblocks; n++) {
    std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[n];
    int schema = blocks_schema_[n];

    if (schema & OPERATOR_DOFS_FACE && schema & OPERATOR_DOFS_CELL) {
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
    } else if (schema == OPERATOR_DOFS_CELL) {
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
* Parallel matvec product A * X.                                              
******************************************************************* */
int Operator::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  // Y = X;
  // return 0;

  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());

  if (X.HasComponent("face")) {
    const Epetra_MultiVector& data = *X.ViewComponent("face");
    for (int f = 0; f < nfaces_owned; f++) Xcopy[f] = data[0][f];
  } 
  if (X.HasComponent("cell")) {
    const Epetra_MultiVector& data = *X.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) Xcopy[c + offset_my_[1]] = data[0][c];
  } 

  preconditioner_->ApplyInverse(Xcopy, Ycopy);

  if (Y.HasComponent("face")) {
    Epetra_MultiVector& data = *Y.ViewComponent("face");
    for (int f = 0; f < nfaces_owned; f++) data[0][f] = Ycopy[f];
  } 
  if (Y.HasComponent("cell")) {
    Epetra_MultiVector& data = *Y.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) data[0][c] = Ycopy[c + offset_my_[1]];
  } 
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



