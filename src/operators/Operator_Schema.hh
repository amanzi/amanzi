/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operator whose unknowns are defined by a general schema.

  The only thing really implemented here is the visitor pattern Op
  acceptors. Everything else should be done in the base class, with
  the exception of special assembly issues.
*/

#ifndef AMANZI_OPERATOR_SCHEMA_HH_
#define AMANZI_OPERATOR_SCHEMA_HH_

#include "DenseVector.hh"

#include "Operator.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Op_Node_Node;
class Op_MeshInjection;

class Operator_Schema : public Operator {
 public:
  // constructors
  // general rectangular operator
  Operator_Schema(const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
                  const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
                  Teuchos::ParameterList& plist,
                  const Schema& schema_row,
                  const Schema& schema_col) :
      Operator(cvs_row, cvs_col, plist, schema_row, schema_col) {
    set_schema_string(schema_col.CreateUniqueName());
  }

  // bijective (square) operator
  Operator_Schema(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                  Teuchos::ParameterList& plist,
                  const Schema& schema) :
      Operator(cvs, cvs, plist, schema, schema) {
    set_schema_string(schema.CreateUniqueName());
  }

  // copy constructor
  virtual Teuchos::RCP<Operator> Clone() const override;

  // required methods
  // -- global methods
  virtual void UpdateRHS(const CompositeVector& source, bool volume_included) override;

  // -- visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Cell_Schema& op,
          const CompositeVector& X, CompositeVector& Y) const override;
  virtual int ApplyMatrixFreeOp(const Op_Face_Schema& op,
          const CompositeVector& X, CompositeVector& Y) const override;
  virtual int ApplyMatrixFreeOp(const Op_Node_Schema& op,
          const CompositeVector& X, CompositeVector& Y) const override;
  virtual int ApplyMatrixFreeOp(const Op_Node_Node& op,
      const CompositeVector& X, CompositeVector& Y) const override;
  virtual int ApplyMatrixFreeOp(const Op_MeshInjection& op,
      const CompositeVector& X, CompositeVector& Y) const override;

  // -- visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;
  virtual void SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;
  virtual void SymbolicAssembleMatrixOp(const Op_Node_Schema& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;
  virtual void SymbolicAssembleMatrixOp(const Op_Node_Node& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;
  virtual void SymbolicAssembleMatrixOp(const Op_MeshInjection& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;

  // -- visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Cell_Schema& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const override;
  virtual void AssembleMatrixOp(const Op_Face_Schema& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const override;
  virtual void AssembleMatrixOp(const Op_Node_Schema& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const override;
  virtual void AssembleMatrixOp(const Op_Node_Node& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const override;
  virtual void AssembleMatrixOp(const Op_MeshInjection& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const override;

  // -- local <-> global communications
  virtual void ExtractVectorCellOp(int c, const Schema& schema,
          WhetStone::DenseVector& v, const CompositeVector& X) const override;
  virtual void AssembleVectorCellOp(int c, const Schema& schema,
          const WhetStone::DenseVector& v, CompositeVector& X) const override;

  virtual void ExtractVectorFaceOp(int f, const Schema& schema,
          WhetStone::DenseVector& v, const CompositeVector& X) const override;
  virtual void AssembleVectorFaceOp(int f, const Schema& schema,
          const WhetStone::DenseVector& v, CompositeVector& X) const override;

  virtual void ExtractVectorNodeOp(int n, const Schema& schema,
          WhetStone::DenseVector& v, const CompositeVector& X) const override;
  virtual void AssembleVectorNodeOp(int n, const Schema& schema,
          const WhetStone::DenseVector& v, CompositeVector& X) const override;

  // debugging methods
  virtual int ApplyAssembled(const CompositeVector& X, CompositeVector& Y, double scalar = 0.0) const override;
};

}  // namespace Operators
}  // namespace Amanzi


#endif



