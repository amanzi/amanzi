/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns/row-col structure are defined by a general schema.

  Supports only the Op_MeshInjection visitor overloads, for a genuinely
  rectangular, cross-mesh, one-nonzero-per-row coupling block (e.g. the
  off-diagonal blocks of a 2x2 TreeOperator coupling two PKs defined on
  different meshes).

  The only thing implemented here is the visitor pattern Op acceptors.
  Everything else is done in the base Operator class.
*/

#ifndef AMANZI_OPERATOR_SCHEMA_HH_
#define AMANZI_OPERATOR_SCHEMA_HH_

#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Operator_Schema : public Operator {
 public:
  // general rectangular operator (domain may differ from range)
  Operator_Schema(const Teuchos::RCP<const CompositeSpace>& cvs_row,
                  const Teuchos::RCP<const CompositeSpace>& cvs_col,
                  Teuchos::ParameterList& plist,
                  const Schema& schema_row,
                  const Schema& schema_col)
    : Operator(cvs_row, cvs_col, plist, schema_row, schema_col)
  {
    set_schema_string("SCHEMA");
  }

  // bijective operator (domain = range)
  Operator_Schema(const Teuchos::RCP<const CompositeSpace>& cvs,
                  Teuchos::ParameterList& plist,
                  const Schema& schema)
    : Operator(cvs, plist, schema)
  {
    set_schema_string("SCHEMA");
  }

  // visit methods for Apply
  virtual int
  ApplyMatrixFreeOp(const Op_MeshInjection& op, const CompositeVector& X, CompositeVector& Y) const override;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_MeshInjection& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const override;

  // visit methods for assemble
  virtual void AssembleMatrixOp(const Op_MeshInjection& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const override;
};

} // namespace Operators
} // namespace Amanzi

#endif
