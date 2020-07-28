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

#ifndef AMANZI_OPERATOR_CELLFACE_SFF_HH_
#define AMANZI_OPERATOR_CELLFACE_SFF_HH_

#include "DenseMatrix.hh"
#include "Operator_FaceCell.hh"
#include "InverseSchurComplement.hh"

namespace Amanzi {
namespace Operators {

class Operator_FaceCellSff : public Operator_FaceCell {
 public:
  // main constructor
  //   The CVS is the domain and range of the operator
  Operator_FaceCellSff(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                       Teuchos::ParameterList& plist) :
      Operator_FaceCell(cvs, plist) {
    // changing schema for the Schur complement
    int schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE;
    schema_col_.Init(schema);
    schema_row_.Init(schema);
    set_schema_string("FACE+CELL Schur to FACE");
  }

  virtual void InitializeInverse(Teuchos::ParameterList& solver_list) override;
  
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

  // visit method for sparsity structure of Schur complement
  // handled in Schur complement -- no cell dofs.
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const override {};

 protected:
  mutable std::vector<Teuchos::RCP<Op_Cell_Face> > schur_ops_;
  Teuchos::RCP<AmanziSolvers::InverseSchurComplement> schur_inv_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

    

