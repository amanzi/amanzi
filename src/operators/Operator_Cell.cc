/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are CELLs.
*/

#include "AmanziTypes.hh"
#include "DenseMatrix.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Face_Cell.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Operator_Cell.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply a source which may or may not have cell volume included already. 
****************************************************************** */
void Operator_Cell::UpdateRHS(const CompositeVector& source,
                              bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    auto rhs_c = rhs_->ViewComponent("cell", false);
    const auto source_c = source.ViewComponent("cell", false);
    for (int c = 0; c != ncells_owned; ++c) {
      rhs_c(0,c) += source_c(0,c) * mesh_->cell_volume(c);
    }
  }
}


/* ******************************************************************
* Visit methods for Apply.
* Apply the local matrices directly as schema is a subset of 
* assembled schema.
****************************************************************** */
int Operator_Cell::ApplyMatrixFreeOp(const Op_Cell_Cell& op,
                                     const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.diag->getLocalLength() == ncells_owned);
  auto Xc = X.ViewComponent("cell");
  auto Yc = Y.ViewComponent("cell");
  const auto dv = op.diag->getLocalViewDevice(); 

  Kokkos::parallel_for(
     "Operator_Cell::ApplyMatrixFreeOp Op_Cell_Cell",
      Xc.extent(0),
      KOKKOS_LAMBDA(const int c) {
        Yc(c,0) += Xc(c,0) * dv(c,0);
      });

  return 0;
}

/* ******************************************************************
* Apply the local matrices directly as schema is a subset of
* assembled schema
****************************************************************** */
int Operator_Cell::ApplyMatrixFreeOp(const Op_Face_Cell& op,
                                     const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.A.size() == nfaces_owned);
  auto Yc = Y.ViewComponent("cell", true);
  auto Xc = X.ViewComponent("cell", true);

  const auto& mc = mesh_->getMeshCache(); 
  
  // Allocate the first time 
  if (op.v.size() != op.A.size()) {
    op.PreallocateWorkVectors();
  }

  auto local_A = op.A; 
  auto local_Av = op.Av; 
  auto local_v = op.v; 

  // parallel matrix-vector product
  Kokkos::parallel_for(
      "Operator_Cell::ApplyMatrixFreeOp Op_Face_Cell COMPUTE",
      nfaces_owned,
      KOKKOS_LAMBDA(const int f) {
        AmanziMesh::Entity_ID_View cells;
        mc.face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);

        int ncells = cells.extent(0);
        auto lv = local_v[f];
        auto lAv = local_Av[f]; 
        auto lA = local_A[f];  
        for (int n = 0; n != ncells; ++n) {
          lv(n) = Xc(cells[n],0);
        }
        lA.Multiply(lv,lAv,false);

        for (int n = 0; n != ncells; ++n) {
          Kokkos::atomic_add(&Yc(cells[n],0), lAv(n));
        }
      });

  return 0;
}

/* ******************************************************************
* Visit methods for symbolic assemble.
* Insert cell-based diagonal matrix.
****************************************************************** */
void Operator_Cell::SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  const auto cell_row_inds = map.GhostIndices<MirrorHost>(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices<MirrorHost>(my_block_col, "cell", 0);

  for (int c = 0; c != ncells_owned; ++c) {
    int row = cell_row_inds[c];
    int col = cell_col_inds[c];

    graph.insertLocalIndices(row, 1, &col);
  }
}


/* ******************************************************************
* Insert each face neighboring cells (ELEMENT/BASE=face, DOFs=cell)
****************************************************************** */
void Operator_Cell::SymbolicAssembleMatrixOp(const Op_Face_Cell& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  std::cout<<"Operator_Cell:SymbolicAssembleMatrixOp"<<std::endl;
  AMANZI_ASSERT(op.A.size() == nfaces_owned);

  int lid_r[2];
  int lid_c[2];

  const auto cell_row_inds = map.GhostIndices<MirrorHost>(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices<MirrorHost>(my_block_col, "cell", 0);

  Kokkos::View<AmanziMesh::Entity_ID*,Kokkos::HostSpace> cells; 
  for (int f = 0; f != nfaces_owned; ++f) {
    op.mesh->face_get_cells_host(f, AmanziMesh::Parallel_type::ALL, cells);
    
    int ncells = cells.size();
    for (int n = 0; n != ncells; ++n) {
      lid_r[n] = cell_row_inds[cells[n]];
      lid_c[n] = cell_col_inds[cells[n]];
    }

    graph.insertLocalIndices(ncells, lid_r, ncells, lid_c);
  }
}


/* ******************************************************************
* Visit methods for assemble
* Insert each cells neighboring cells.
****************************************************************** */
void Operator_Cell::AssembleMatrixOp(const Op_Cell_Cell& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  std::cout<<"Operator_Cell::AssembleMatrixOp OP_CELL_CELL"<<std::endl;
  AMANZI_ASSERT(op.diag->getNumVectors() == 1);

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);
  const auto dv = op.diag->getLocalViewDevice(); 

  // hard-coded version, interfaces TBD...
  auto proc_mat = mat.getLocalMatrix();
  auto offproc_mat = mat.getOffProcLocalMatrix();
  int nrows_local = mat.getMatrix()->getNodeNumRows();
  
  Kokkos::parallel_for(
      "Operator_Cell::AssembleMatrixOp::Cell_Cell",
      ncells_owned,
      KOKKOS_LAMBDA(const int& c) {
        if (cell_row_inds(c) < nrows_local) {
          proc_mat.sumIntoValues(cell_row_inds(c),
                  &cell_col_inds(c), 1, &dv(c,0), true, false);
        } else {
          offproc_mat.sumIntoValues(cell_row_inds(c) - nrows_local,
                  &cell_col_inds(c), 1, &dv(c,0), true, false);
        }
      });
}


void Operator_Cell::AssembleMatrixOp(const Op_Face_Cell& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.A.size() == nfaces_owned);

  op.A.update_entries_host(); 

  #ifdef OUTPUT_CUDA 
  std::cout<<std::endl;
  std::cout<<"host: "; 
  for(int i = 0 ; i < nfaces_owned; ++i){
    std::cout<<op.A.at_host(i)(0,0)<<" - ";
  }
  std::cout<<std::endl;

  //op.A.update_entries_device(); 

  std::cout<<"device: ";
  Kokkos::parallel_for(
    "", 
    nfaces_owned, 
    KOKKOS_LAMBDA(const int i){ 
      printf("%.4f - ",op.A[i](0,0)); 
    });
    Kokkos::fence(); 
  std::cout<<std::endl;
  std::cout<<"Operator_Cell::AssembleMatrixOp"<<std::endl;
  #endif 
  const auto cell_row_inds = map.GhostIndices<>(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices<>(my_block_col, "cell", 0);

  // hard-coded version, interfaces TBD...
  auto proc_mat = mat.getLocalMatrix();
  auto offproc_mat = mat.getOffProcLocalMatrix();
  int nrows_local = mat.getMatrix()->getNodeNumRows();

  const AmanziMesh::Mesh* mesh = op.mesh.get();
  #ifdef OUTPUT_CUDA 
  std::cout<<"nfaces_owned: "<<nfaces_owned<<std::endl;
  #endif 
  Kokkos::parallel_for(
      "Operator_Cell::AssembleMatrixOp::Face_Cell",
      nfaces_owned,
      KOKKOS_LAMBDA(const int f) {
        AmanziMesh::Entity_ID_View cells;
        mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
    
        auto A_f = op.A[f];
        int nc = cells.extent(0); 
        for (int n = 0; n < nc; ++n) {
          if (cell_row_inds[cells[n]] < nrows_local) {
            for (int m = 0; m < cells.extent(0); ++m) {
              #ifdef OUTPUT_CUDA 
              printf("# %d = %.4f %lld ( %d / %d / %d / %d) %d,%d  - ",
                         f,A_f(n,m),cells.extent(0),n,m,cells(n),
                         cells(m),cell_row_inds(cells(n)),cell_col_inds(cells(m))); 
              #endif 
              proc_mat.sumIntoValues(cell_row_inds(cells(n)),
                      &cell_col_inds(cells(m)), 1,
                      &A_f(n,m), true, true);
            }
          } else {
            for (int m = 0; m != cells.extent(0); ++m) {
              offproc_mat.sumIntoValues(cell_row_inds(cells(n))-nrows_local,
                      &cell_col_inds(cells(m)), 1,
                      &A_f(n,m), true, true);
            }
          }
        }
      });
  #ifdef OUTPUT_CUDA 
  Kokkos::fence(); 
  std::cout<<std::endl;

  std::cout<<"device after: ";
  Kokkos::parallel_for(
    "", 
    nfaces_owned, 
    KOKKOS_LAMBDA(const int i){ 
      printf("(%d) %.4f - ",i,op.A[i](0,0)); 
    });
    Kokkos::fence(); 
  std::cout<<std::endl;
  #endif 

}


}  // namespace Operators
}  // namespace Amanzi
