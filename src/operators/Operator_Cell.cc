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

  AMANZI_ASSERT(op.csr.size() == nfaces_owned);
  auto Yc = Y.ViewComponent("cell", true);
  auto Xc = X.ViewComponent("cell", true);

  CSR_Matrix local_csr = op.csr; 

  const AmanziMesh::Mesh* mesh = mesh_.get();

  // Allocate for first time 
  if(op.csr_v_.size() != local_csr.size()){

    op.csr_v_ = CSR<double,1,DeviceOnlyMemorySpace>(local_csr.size());
    op.csr_Av_ = CSR<double,1,DeviceOnlyMemorySpace>(local_csr.size()); 

    int total1 = 0; 
    int total2 = 0; 
    // CSR version 
    // 1. Compute size 
    for (int i=0; i!=local_csr.size(); ++i) {
      total1 += local_csr.size(i,0);
      total2 += local_csr.size(i,1);
    }

    Kokkos::parallel_for(
      "Operator_Cell::ApplyMatrixFreeOp Op_Face_Cell COPY",
      local_csr.size(),
      KOKKOS_LAMBDA(const int& i){      
        op.csr_v_.sizes_(i,0) = local_csr.size(i,0);
        op.csr_v_.row_map_(i) = local_csr.size(i,0);
        op.csr_Av_.sizes_(i,0) = local_csr.size(i,1);
        op.csr_Av_.row_map_(i) = local_csr.size(i,1);
      });

    op.csr_v_.prefix_sum_device(total1); 
    op.csr_Av_.prefix_sum_device(total2); 
  }

  CSR<double,1,DeviceOnlyMemorySpace> csr_v = op.csr_v_;
  CSR<double,1,DeviceOnlyMemorySpace> csr_Av = op.csr_Av_;

  Kokkos::parallel_for(
      "Operator_Cell::ApplyMatrixFreeOp Op_Face_Cell COMPUTE",
      op.csr.size(),
      KOKKOS_LAMBDA(const int f) {
        AmanziMesh::Entity_ID_View cells;
        mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);

        int ncells = cells.extent(0);

        WhetStone::DenseVector<DeviceOnlyMemorySpace> vv(
          csr_v.at(f),csr_v.size(f));
        WhetStone::DenseVector<DeviceOnlyMemorySpace> Avv(
          csr_Av.at(f), csr_Av.size(f));

        for (int n = 0; n != ncells; ++n) {
          vv(n) = Xc(cells[n],0);
        }
        WhetStone::DenseMatrix lm(
          local_csr.at(f),
          local_csr.size(f,0),local_csr.size(f,1)); 
        lm.Multiply(vv,Avv,false);

        for (int n = 0; n != ncells; ++n) {
          Kokkos::atomic_add(&Yc(cells[n],0), Avv(n));

          // if (cells[n] == 0 || cells[n] == 9) {
          //   std::cout << std::setprecision(16) << "Apply at f(" << f << "): v = " << vv(0);
          //   if (csr_Av.size(f) > 1) std::cout << std::setprecision(16)  << "," << vv(1);
          //   std::cout << std::setprecision(16)  << ", A = " << lm << ", Av = " << Avv(0);
          //   if (csr_Av.size(f) > 1) std::cout << std::setprecision(16)  << "," << Avv(1);
          //   std::cout << std::endl;
          // }

        }
      });

  return 0;
}

#if 0 
/* ******************************************************************
* Visit methods for symbolic assemble.
* Insert cell-based diagonal matrix.
****************************************************************** */
void Operator_Cell::SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    int row = cell_row_inds[c];
    int col = cell_col_inds[c];

    ierr |= graph.InsertMyIndices(row, 1, &col);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Insert each face neighboring cells (ELEMENT/BASE=face, DOFs=cell)
****************************************************************** */
void Operator_Cell::SymbolicAssembleMatrixOp(const Op_Face_Cell& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  std::vector<int> lid_r;
  std::vector<int> lid_c;

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    AmanziMesh::Entity_ID_View cells; 
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
    
    int ncells = cells.size();
    lid_r.resize(ncells);
    lid_c.resize(ncells);    
    for (int n = 0; n != ncells; ++n) {
      lid_r[n] = cell_row_inds[cells[n]];
      lid_c[n] = cell_col_inds[cells[n]];
    }

    ierr |= graph.InsertMyIndices(ncells, lid_r.data(), ncells, lid_c.data());
  }

  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert each cells neighboring cells.
****************************************************************** */
void Operator_Cell::AssembleMatrixOp(const Op_Cell_Cell& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.diag.getNumVectors() == 1);

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);
  const auto dv = op.diag->getLocalViewDevice(); 

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    int row = cell_row_inds[c];
    int col = cell_col_inds[c];

    ierr |= mat.SumIntoMyValues(row, 1, &dv(0,c), &col);
  }
  AMANZI_ASSERT(!ierr);
}


void Operator_Cell::AssembleMatrixOp(const Op_Face_Cell& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);
  
  // ELEMENT: face, DOF: cell
  std::vector<int> lid_r;
  std::vector<int> lid_c;

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    AmanziMesh::Entity_ID_View cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
    
    int ncells = cells.size();
    lid_r.resize(ncells);
    lid_c.resize(ncells);    
    for (int n = 0; n != ncells; ++n) {      
      lid_r[n] = cell_row_inds[cells[n]];
      lid_c[n] = cell_col_inds[cells[n]];
    }
   
    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[f]);
    AMANZI_ASSERT(ierr>=0);
  }
  AMANZI_ASSERT(ierr>=0);
}
#endif 

}  // namespace Operators
}  // namespace Amanzi
