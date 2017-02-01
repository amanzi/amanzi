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

#include "Op_Node_Node.hh"
#include "Operator.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Operator_Schema : public Operator {
 public:
  // main constructor
  // The input CVS is the domain and range of the operator.
  Operator_Schema(const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
                  const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
                  Teuchos::ParameterList& plist,
                  const Schema& schema_row,
                  const Schema& schema_col) :
      Operator(cvs_row, cvs_col, plist, schema_row, schema_col) {
    set_schema_string(schema_col.CreateUniqueName());
  }

  // required methods
  // -- global methods
  virtual void SymbolicAssembleMatrix();
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;
  virtual void UpdateRHS(const CompositeVector& source, bool volume_included);

  // -- visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Cell_Schema& op,
          const CompositeVector& X, CompositeVector& Y) const;
  virtual int ApplyMatrixFreeOp(const Op_Node_Node& op,
      const CompositeVector& X, CompositeVector& Y) const override;

  // -- visit methods for ApplyTranspose 
  virtual int ApplyTransposeMatrixFreeOp(const Op_Cell_Schema& op,
          const CompositeVector& X, CompositeVector& Y) const;

  // -- visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const;
  virtual void SymbolicAssembleMatrixOp(const Op_Node_Node& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;
  
  // -- visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Cell_Schema& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const;
  virtual void AssembleMatrixOp(const Op_Node_Node& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const;

  // -- local <-> global communications
  virtual void ExtractVectorOp(int c, const Schema& schema,
          WhetStone::DenseVector& v, const CompositeVector& X) const;
  virtual void AssembleVectorOp(int c, const Schema& schema,
          const WhetStone::DenseVector& v, CompositeVector& X) const;

  // debugging methods
  int ApplyAssembled(const CompositeVector& X, CompositeVector& Y, double scalar = 0.0) const;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

    

