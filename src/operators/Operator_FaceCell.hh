/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Operator whose unknowns are CELL + FACE

  The only thing really implemented here is the visitor pattern Op
  acceptors. Everything else should be done in the base class, with
  the exception of special assembly issues.
*/

#ifndef AMANZI_OPERATOR_CELLFACE_HH_
#define AMANZI_OPERATOR_CELLFACE_HH_

#include "DenseMatrix.hh"
#include "Operator_Cell.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class Operator_FaceCell : public Operator_Cell {
 public:
  // main constructor
  // The input CVS is the domain and range of the operator.
  Operator_FaceCell(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                    Teuchos::ParameterList& plist)
    : Operator_Cell(cvs, plist, OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)
  {
    set_schema_string("FACE+CELL");
  }

  virtual Teuchos::RCP<Operator> Clone() const;

  // visit methods for Apply
  using Operator_Cell::ApplyMatrixFreeOp;
  virtual int
  ApplyMatrixFreeOp(const Op_Cell_FaceCell& op, const CompositeVector& X, CompositeVector& Y) const;
  virtual int
  ApplyMatrixFreeOp(const Op_Cell_Face& op, const CompositeVector& X, CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const;

  // visit methods for symbolic assemble
  using Operator_Cell::SymbolicAssembleMatrixOp;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  // visit methods for actual assemble
  using Operator_Cell::AssembleMatrixOp;
  virtual void AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Cell_Face& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;
};

} // namespace Operators
} // namespace Amanzi

#endif
