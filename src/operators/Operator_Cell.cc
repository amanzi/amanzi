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

  CSR<double,1,DeviceOnlyMemorySpace> csr_v = op.csr_v_;
  CSR<double,1,DeviceOnlyMemorySpace> csr_Av = op.csr_Av_;

  // Allocate for first time 
  if(csr_v.size_host() != local_csr.size_host()){

    csr_v.set_size(local_csr.size_host());
    csr_Av.set_size(local_csr.size_host());  

    int total1 = 0; 
    int total2 = 0; 
    // CSR version 
    // 1. Compute size 
    for (int i=0; i!=local_csr.size_host(); ++i) {
      total1 += local_csr.size_host(i,0);
      total2 += local_csr.size_host(i,1);
    }

    Kokkos::parallel_for(
      "Operator_Cell::ApplyMatrixFreeOp Op_Face_Cell COPY",
      local_csr.size_host(),
      KOKKOS_LAMBDA(const int& i){      
        csr_v.sizes_.view_device()(i,0) = local_csr.size(i,0);
        csr_v.row_map_.view_device()(i) = local_csr.size(i,0);
        csr_Av.sizes_.view_device()(i,0) = local_csr.size(i,1);
        csr_Av.row_map_.view_device()(i) = local_csr.size(i,1);
      });

    csr_v.prefix_sum_device(total1); 
    csr_Av.prefix_sum_device(total2); 
  }

  Kokkos::parallel_for(
      "Operator_Cell::ApplyMatrixFreeOp Op_Face_Cell COMPUTE",
      op.csr.size_host(),
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
        WhetStone::DenseMatrix<DeviceOnlyMemorySpace> lm(
          local_csr.at(f),
          local_csr.size(f,0),local_csr.size(f,1)); 
        lm.Multiply(vv,Avv,false);

        for (int n = 0; n != ncells; ++n) {
          Kokkos::atomic_add(&Yc(cells[n],0), Avv(n));
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
  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

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
  AMANZI_ASSERT(op.csr.size() == nfaces_owned);

  std::vector<int> lid_r;
  std::vector<int> lid_c;

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

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

    graph.insertLocalIndices(ncells, lid_r.data(), ncells, lid_c.data());
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
  AMANZI_ASSERT(op.diag->getNumVectors() == 1);

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);
  const auto dv = op.diag->getLocalViewDevice(); 

  for (int c = 0; c != ncells_owned; ++c) {
    int row = cell_row_inds[c];
    int col = cell_col_inds[c];
    mat.sumIntoLocalValues(row, 1, &dv(0,c), &col);
  }
}


void Operator_Cell::AssembleMatrixOp(const Op_Face_Cell& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.csr.size() == nfaces_owned);
  
  // ELEMENT: face, DOF: cell
  std::vector<int> lid_r;
  std::vector<int> lid_c;

  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);
  // need to make this const
  // op.csr.update_entries_host();

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

    // what we would like to do...
    //WhetStone::DenseMatrix<Amanzi::DeviceOnlyMemorySpace> m(op.csr.at(f), op.csr.size(f,0), op.csr.size(f,1));
    //mat.sumIntoLocalValues(lid_r.data(), lid_c.data(), m_host);

    // what we can do
    WhetStone::DenseMatrix<Amanzi::DefaultHostMemorySpace> lm(
            op.csr.at_host(f),
            op.csr.size_host(f,0),op.csr.size_host(f,1)); 
    mat.sumIntoLocalValues((const int*) lid_r.data(), (const int*) lid_c.data(), lm);
    std::cout << "OpCell add: " << lm(0,0) << std::endl;
  }
}


}  // namespace Operators
}  // namespace Amanzi
