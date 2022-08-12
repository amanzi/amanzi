/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky(dasvyat@lanl.gov)

  Operator whose unknowns are CELLs and BOUNDARY FACES

  The only thing really implemented here is the visitor pattern Op
  acceptors. Everything else should be done in the base class, with
  the exception of special assembly issues.
*/

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
  Operator_CellBndFace(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                Teuchos::ParameterList& plist,
                int schema) :
    Operator_Cell(cvs, plist, schema) {
    set_schema_string("CELLBNDFACE");
  }

  // visit methods for apply
  using Operator_Cell::ApplyMatrixFreeOp;
  virtual int ApplyMatrixFreeOp(const Op_Face_CellBndFace& op,
                                const CompositeVector& X, CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                                const CompositeVector& X, CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                const CompositeVector& X, CompositeVector& Y) const;

  // visit methods for symbolic assemble
  using Operator_Cell::SymbolicAssembleMatrixOp;
  virtual void SymbolicAssembleMatrixOp(const Op_Face_CellBndFace& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const;
  
  virtual void SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                        const SuperMap& map, GraphFE& graph,
                                        int my_block_row, int my_block_col) const;
  
  // visit methods for assemble
  using Operator_Cell::AssembleMatrixOp;
  virtual void AssembleMatrixOp(const Op_Face_CellBndFace& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                const SuperMap& map, MatrixFE& mat,
                                int my_block_row, int my_block_col) const;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

    

