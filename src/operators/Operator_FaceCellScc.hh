/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_CELLFACE_SCC_HH_
#define AMANZI_OPERATOR_CELLFACE_SCC_HH_

#include "DenseMatrix.hh"
#include "Operator_Cell.hh"

/* ******************************************************************
Operator whose unknowns are CELL + FACE, but which assembles the CELL only
system and Schur complements the face.

This uses special assembly.  Apply is done as if we had the full FACE+CELL
system.  SymbolicAssembly() is done as if we had the CELL system, but with an
additional step to get the layout due to the Schur'd system on FACE+CELL.
Assemble, however, is done using a totally different approach.

---------------------------------------------------------------------

1. Operator is a linear operator acting from linear space X to linear
space Y. These spaces are described by CompositeVectors (CV). A few
maps X->Y is supported. 

At the moment X = Y. Extension to TreeVectors should not be done in 
this class.

2. Operator is an un-ordered additive collection of lower-rank (or 
equal) simple operators. During its construction, an operator can 
only grow by assimilating more operators. 

At the moment, an operator cannot be split into two operators, but
there are no desing restriction for doing it in the future.

3. A simple operator (a set of 1 operators) is defined by triple:
scheme + elemental matrices + diagonal matrix. The schema specifies
structure of elemental matrices, e.g. cell-based matrices 
representing interation between face-based unknowns.

4. Operator can be converted to Epetra_FECrsMatrix matrix to generate
a preconditioner. This operation cannot be applied to a subset of
defining operators. 
 
Note. The operators can be initialized from other operators.
    Since data are never copied by default, we have to track 
    down the ownership of data.
****************************************************************** */ 

namespace Amanzi {
namespace Operators {

class Operator_FaceCellScc : public Operator_Cell {
 public:
  // main constructor
  // The input CVS is the domain and range of the operator
  Operator_FaceCellScc(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                       Teuchos::ParameterList& plist) :
      Operator_Cell(cvs, plist, OPERATOR_SCHEMA_DOFS_CELL) {
    set_schema_string("FACE+CELL Schur to CELL");
  }

  // Special Apply Inverse required to deal with schur complement
  virtual int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;

  // Special AssembleMatrix required to deal with schur complement
  virtual void AssembleMatrix(const SuperMap& map,
          MatrixFE& matrix, int my_block_row, int my_block_col) const;
  
  // visit method for Apply -- this is identical to Operator_FaceCell's
  // version.
  virtual int ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
      const CompositeVector& X, CompositeVector& Y) const;

  // driver symbolic assemble creates the face-only supermap
  virtual void SymbolicAssembleMatrix();

  // visit method for sparsity structure of Schur complement
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const;

 protected:
  mutable std::vector<Teuchos::RCP<Op_Cell_Cell> > diag_ops_;
  mutable std::vector<Teuchos::RCP<Op_Face_Cell> > schur_ops_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

    

