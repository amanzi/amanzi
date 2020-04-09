/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_NODE_HH_
#define AMANZI_OPERATOR_NODE_HH_

#include "DenseMatrix.hh"
#include "OperatorDefs.hh"
#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Operator_Node : public Operator {
 public:
  // main constructor
  //   The CVS is the domain and range of the operator
  Operator_Node(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                Teuchos::ParameterList& plist)
    : Operator(cvs, plist, OPERATOR_SCHEMA_DOFS_NODE)
  {
    set_schema_string("NODE");
    cell_max_nodes = mesh_->cell_get_max_nodes();
  }

  // rhs update which multiplies by cell
  virtual void
  UpdateRHS(const CompositeVector& source, bool volume_included) override;

  // visit methods for Apply
  virtual int
  ApplyMatrixFreeOp(const Op_Cell_Node& op, const CompositeVector& X,
                    CompositeVector& Y) const override;

  virtual int
  ApplyMatrixFreeOp(const Op_Node_Node& op, const CompositeVector& X,
                    CompositeVector& Y) const override;

  // visit methods for symbolic assemble
  virtual void
  SymbolicAssembleMatrixOp(const Op_Cell_Node& op, const SuperMap& map,
                           GraphFE& graph, int my_block_row,
                           int my_block_col) const override;

  virtual void
  SymbolicAssembleMatrixOp(const Op_Node_Node& op, const SuperMap& map,
                           GraphFE& graph, int my_block_row,
                           int my_block_col) const override;

  // visit methods for assemble
  virtual void
  AssembleMatrixOp(const Op_Cell_Node& op, const SuperMap& map, MatrixFE& mat,
                   int my_block_row, int my_block_col) const override;

  virtual void
  AssembleMatrixOp(const Op_Node_Node& op, const SuperMap& map, MatrixFE& mat,
                   int my_block_row, int my_block_col) const override;

 protected:
  int cell_max_nodes;
};

} // namespace Operators
} // namespace Amanzi

#endif
