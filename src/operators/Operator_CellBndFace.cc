/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky(dasvyat@lanl.gov) 

  Operator whose unknowns are CELLs and BOUNDARY FACES
*/

#include "DenseMatrix.hh"
//#include "Op_Cell_Cell.hh"
#include "Op_Face_CellFace.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Operator_CellBndFace.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply a source which may or may not have cell volume included already. 
****************************************************************** */
void Operator_CellBndFace::UpdateRHS(const CompositeVector& source,
                              bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    Epetra_MultiVector& rhs_c = *rhs_->ViewComponent("cell", false);
    const Epetra_MultiVector& source_c = *source.ViewComponent("cell", false);
    for (int c = 0; c != ncells_owned; ++c) {
      rhs_c[0][c] += source_c[0][c] * mesh_->cell_volume(c);
    }
  }
}


/* ******************************************************************
* Apply the local matrices directly as schema is a subset of
* assembled schema
****************************************************************** */
int Operator_CellBndFace::ApplyMatrixFreeOp(const Op_Face_CellFace& op,
                                     const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == nfaces_owned);
  
  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);
  const Epetra_MultiVector& Xbnd = *X.ViewComponent("boundary face", true);

  Y.PutScalarGhosted(0.);
  Epetra_MultiVector& Yc = *Y.ViewComponent("cell", true);
  Epetra_MultiVector& Ybnd = *Y.ViewComponent("boundary face", true); 

  AmanziMesh::Entity_ID_List cells;
  for (int f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 2){
      WhetStone::DenseVector v(ncells), av(ncells);
      for (int n=0; n!=ncells; ++n) {
        v(n) = Xc[0][cells[n]];
      }

      const WhetStone::DenseMatrix& Aface = op.matrices[f];
      Aface.Multiply(v, av, false);

      for (int n=0; n!=ncells; ++n) {
        Yc[0][cells[n]] += av(n);
      }
    }else if (ncells==1){
      int bf = mesh_->exterior_face_map(false).LID(mesh_->face_map(false).GID(f));

      WhetStone::DenseVector v(2), av(2);
      v(0) = Xc[0][cells[0]];
      v(1) = Xbnd[0][bf];
      
      const WhetStone::DenseMatrix& Aface = op.matrices[f];
      Aface.Multiply(v, av, false);

      Yc[0][cells[0]] += av(0);
      Ybnd[0][bf] += av(1);    
    }
  }

  Y.GatherGhostedToMaster("cell",Add);
  Y.GatherGhostedToMaster("boundary face",Add);
  
  return 0;
}



/* ******************************************************************
* Insert each cells neighboring cells.
****************************************************************** */
void Operator_CellBndFace::SymbolicAssembleMatrixOp(const Op_Face_CellFace& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  // ELEMENT: face, DOF: cell, bnd_face
  int lid_r[2];
  int lid_c[2];
  const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
  const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);
  const std::vector<int>& bndface_row_inds = map.GhostIndices("boundary face", my_block_row);
  const std::vector<int>& bndface_col_inds = map.GhostIndices("boundary face", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List cells;
  for (int f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    
    int ncells = cells.size();
    if (ncells == 2){
      for (int n=0; n!=ncells; ++n) {
        lid_r[n] = cell_row_inds[cells[n]];
        lid_c[n] = cell_col_inds[cells[n]];
      }
    }else if (ncells==1){
      lid_r[0] = cell_row_inds[cells[0]];
      lid_c[0] = cell_col_inds[cells[0]];
      int bf = mesh_->exterior_face_map(false).LID(mesh_->face_map(false).GID(f));
      lid_r[1] = bndface_row_inds[bf];
      lid_c[1] = bndface_col_inds[bf];
    }


    ierr |= graph.InsertMyIndices(2, lid_r, 2, lid_c);
  }
  ASSERT(!ierr);
}



void Operator_CellBndFace::AssembleMatrixOp(const Op_Face_CellFace& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == nfaces_owned);
  
  // ELEMENT: face, DOF: cell,  bnd_face
  int lid_r[2];
  int lid_c[2];
  const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
  const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);
  const std::vector<int>& bndface_row_inds = map.GhostIndices("boundary face", my_block_row);
  const std::vector<int>& bndface_col_inds = map.GhostIndices("boundary face", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List cells;
  for (int f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    
    int ncells = cells.size();
    if (ncells == 2){
      for (int n=0; n!=ncells; ++n) {
        lid_r[n] = cell_row_inds[cells[n]];
        lid_c[n] = cell_col_inds[cells[n]];
      }
    }else if (ncells==1){
      lid_r[0] = cell_row_inds[cells[0]];
      lid_c[0] = cell_col_inds[cells[0]];
      int bf = mesh_->exterior_face_map(false).LID(mesh_->face_map(false).GID(f));
      lid_r[1] = bndface_row_inds[bf];
      lid_c[1] = bndface_col_inds[bf];
      // std::cout  << "bnd "<<lid_r[0]<<" "<<lid_r[1]<<"\n";
      //std::cout<<op.matrices[f]<<"\n";
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[f]);
    ASSERT(!ierr);
  }
  ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi

