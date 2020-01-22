/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Daniil Svyatsky(dasvyat@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_WITH_CELLBND_HH_
#define AMANZI_OPERATOR_WITH_CELLBND_HH_

#include "DenseMatrix.hh"
#include "Operator_Cell.hh"

namespace Amanzi {
namespace Operators {

class Operator_CellBndFace : public Operator_Cell {
 public:
  // main constructor
  //   The CVS is the domain and range of the operator
  Operator_CellBndFace(const Teuchos::RCP<const CompositeSpace>& cvs,
                       Teuchos::ParameterList& plist, int schema)
    : Operator_Cell(cvs, plist, schema)
  {
    set_schema_string("CELLBNDFACE");
  }

  // visit methods for apply
  virtual int
  ApplyMatrixFreeOp(const Op_Face_CellBndFace& op, const CompositeVector& X,
                    CompositeVector& Y) const;
  
  virtual void getLocalDiagCopy(CompositeVector& X) const;

  // virtual int
  // ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
  //                   const CompositeVector& X, CompositeVector& Y) const;
  // virtual int
  // ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
  //                   const CompositeVector& X, CompositeVector& Y) const;

  // // visit methods for symbolic assemble
  // virtual void
  // SymbolicAssembleMatrixOp(const Op_Face_CellBndFace& op, const SuperMap& map,
  //                          GraphFE& graph, int my_block_row,
  //                          int my_block_col) const;

  // virtual void
  // SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
  //                          const SuperMap& map, GraphFE& graph,
  //                          int my_block_row, int my_block_col) const;
  // virtual void
  // SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
  //                          const SuperMap& map, GraphFE& graph,
  //                          int my_block_row, int my_block_col) const;

  // // visit methods for assemble
  // virtual void
  // AssembleMatrixOp(const Op_Face_CellBndFace& op, const SuperMap& map,
  //                  MatrixFE& mat, int my_block_row, int my_block_col) const;

  // virtual void
  // AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op, const SuperMap& map,
  //                  MatrixFE& mat, int my_block_row, int my_block_col) const;
  // virtual void
  // AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op, const SuperMap& map,
  //                  MatrixFE& mat, int my_block_row, int my_block_col) const;
};

} // namespace Operators
} // namespace Amanzi

#endif
