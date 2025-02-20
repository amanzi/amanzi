/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OP_FACE_FACE_HH_
#define AMANZI_OP_FACE_FACE_HH_

#include <vector>
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Face_Face : public Op {
 public:
  Op_Face_Face(const std::string& name, const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_FACE, name, mesh)
  {
    diag =
      Teuchos::rcp(new Epetra_MultiVector(mesh->getMap(AmanziMesh::Entity_kind::FACE, false), 1));
    diag_shadow =
      Teuchos::rcp(new Epetra_MultiVector(mesh->getMap(AmanziMesh::Entity_kind::FACE, false), 1));
  }

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X, CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const
  {
    assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const
  {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }

  virtual void Rescale(const CompositeVector& scaling)
  {
    if (scaling.HasComponent("face")) {
      const Epetra_MultiVector& s_f = *scaling.ViewComponent("face", false);
      for (int k = 0; k != s_f.NumVectors(); ++k) {
        for (int i = 0; i != s_f.MyLength(); ++i) { (*diag)[k][i] *= s_f[0][i]; }
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi

#endif
