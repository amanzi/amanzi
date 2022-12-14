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

  Operator whose unknowns given by two set of indices of Op which are
  subsets of the domain and range given by composite vector spaces.
  The only thing really implemented here is the visitor pattern Op
  acceptors.
*/

#ifndef AMANZI_OPERATOR_DIAGONAL_HH_
#define AMANZI_OPERATOR_DIAGONAL_HH_

#include "DenseMatrix.hh"
#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Operator_Diagonal : public Operator {
 public:
  // The input CVSs define the domain and range of the operator.
  Operator_Diagonal(const Teuchos::RCP<const CompositeVectorSpace>& cvs_row,
                    const Teuchos::RCP<const CompositeVectorSpace>& cvs_col,
                    Teuchos::ParameterList& plist,
                    int schema)
    : Operator(cvs_row, cvs_col, plist, Schema(schema), Schema(schema))
  {
    row_compname_ = *(cvs_row->begin());
    col_compname_ = *(cvs_col->begin());
    set_schema_string("INDICES");
  }

  // copy constructor
  virtual Teuchos::RCP<Operator> Clone() const override;

  // required methods
  // -- global methods that cannot be aplied to this operator
  // virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const override { AMANZI_ASSERT(false); return 0; }
  // virtual void UpdateRHS(const CompositeVector& source, bool volume_included) override { AMANZI_ASSERT(false); }

  // visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Diagonal& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const override;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Diagonal& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const override;

  // visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Diagonal& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const override;

 private:
  std::string row_compname_, col_compname_;
};

} // namespace Operators
} // namespace Amanzi


#endif
