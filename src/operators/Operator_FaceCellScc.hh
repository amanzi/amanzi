/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are CELL + FACE, but which assembles the
  CELL only system and Schur complements the face.
*/

#ifndef AMANZI_OPERATOR_CELLFACE_SCC_HH_
#define AMANZI_OPERATOR_CELLFACE_SCC_HH_

#include "DenseMatrix.hh"
#include "Operator_Cell.hh"

namespace Amanzi {
namespace Operators {

class Operator_FaceCellScc : public Operator_Cell {
 public:
  // The input CVS is the domain and range of the operator
  Operator_FaceCellScc(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                       Teuchos::ParameterList& plist) :
      Operator_Cell(cvs, plist, OPERATOR_SCHEMA_DOFS_CELL) {
    set_schema_string("FACE+CELL Schur to CELL");
  }

  virtual Teuchos::RCP<Operator> Clone() const override;

  // cannot do an assembled forward apply as the assembled thing is not the full
  // operator
  virtual int Apply(const CompositeVector& X, CompositeVector& Y) const override {
    return Apply(X,Y,0.0);
  }
  virtual int
  Apply(const CompositeVector& X, CompositeVector& Y, double scalar) const override {
    return ApplyUnassembled(X,Y,scalar);
  }
  // Special Apply Inverse required to deal with schur complement
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const override;

  // Special AssembleMatrix required to deal with schur complement
  virtual void AssembleMatrix(const SuperMap& map,
          MatrixFE& matrix, int my_block_row, int my_block_col) const override;

  // visit method for Apply -- this is identical to Operator_FaceCell's
  // version.
  virtual int ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
      const CompositeVector& X, CompositeVector& Y) const override;

  // driver symbolic assemble creates the face-only supermap
  virtual void SymbolicAssembleMatrix() override;

  // visit method for sparsity structure of Schur complement
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override;

 protected:
  mutable std::vector<Teuchos::RCP<Op_Cell_Cell> > diag_ops_;
  mutable std::vector<Teuchos::RCP<Op_Face_Cell> > schur_ops_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

    

