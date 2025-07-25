/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Operator whose unknowns are CELL + FACE but possibly from multiple
  meshes in contact.
*/

#include "DenseMatrix.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
#include "Op_Face_Face.hh"
#include "Op_SurfaceCell_SurfaceCell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
#include "Operator_MultiMesh.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructor
****************************************************************** */
Operator_MultiMesh::Operator_MultiMesh(Teuchos::ParameterList& plist,
                                       Teuchos::RCP<Operator> global_op,
                                       Teuchos::RCP<Op> local_op,
                                       std::vector<int>& interface_block,
                                       InterfaceData& interface_data)
  : Operator(global_op->get_domain_map(), plist, local_op->schema_old()),
    interface_block_(interface_block),
    interface_data_(interface_data)
{
  cell_max_faces = 20; // FIXME
  local_op->schema_string = "Diffusion: FACE_CELL MultiMesh";
  OpPushBack(local_op);
}


/* ******************************************************************
* Visit methods for symbolic assemble: element=cell; dofs={face,cell}
****************************************************************** */
void
Operator_MultiMesh::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                             const SuperMap& map,
                                             GraphFE& graph,
                                             int my_block_row,
                                             int my_block_col) const
{
  std::vector<int> lid_r(2 * cell_max_faces + 1);
  std::vector<int> lid_c(2 * cell_max_faces + 1);

  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);
  const std::vector<int>& cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const std::vector<int>& cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    int k = nfaces + 1;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];

      int other_block = interface_block_[f];
      if (other_block >= 0) {
        const std::vector<int>& other_inds = map.GhostIndices(other_block, "face", 0);
        for (const auto& data : interface_data_[f]) {
          int f2 = data.first;
          lid_r[k] = other_inds[f2];
          lid_c[k] = other_inds[f2];
          k++;
        }
      }
    }

    lid_r[nfaces] = cell_row_inds[c];
    lid_c[nfaces] = cell_col_inds[c];

    ierr |= graph.InsertMyIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: element=cell; dofs={face, cell}.
****************************************************************** */
void
Operator_MultiMesh::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                     const SuperMap& map,
                                     MatrixFE& mat,
                                     int my_block_row,
                                     int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(2 * cell_max_faces + 1);
  std::vector<int> lid_c(2 * cell_max_faces + 1);

  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);
  const std::vector<int>& cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const std::vector<int>& cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    int k = nfaces + 1;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];

      int other_block = interface_block_[f];
      if (other_block >= 0) {
        const std::vector<int>& other_inds = map.GhostIndices(other_block, "face", 0);
        for (const auto& data : interface_data_[f]) {
          int f2 = data.first;
          lid_r[k] = other_inds[f2];
          lid_c[k] = other_inds[f2];
          k++;
        }
      }
    }

    lid_r[nfaces] = cell_row_inds[c];
    lid_c[nfaces] = cell_col_inds[c];

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator>
Operator_MultiMesh::Clone() const
{
  return Teuchos::rcp(new Operator_MultiMesh(*this));
}

} // namespace Operators
} // namespace Amanzi
