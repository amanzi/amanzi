/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_EDGE_HH_
#define AMANZI_OPERATOR_EDGE_HH_

#include "DenseMatrix.hh"
#include "Operator.hh"

/* ******************************************************************
Operator whose unknowns are EDGE

NOTE that the only thing really implemented here is the visitor pattern Op
acceptors.  Everything else should be done in the base class, with the
exception of special assembly issues.

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

class Operator_Edge : public Operator {
 public:
  // main constructor
  // The input CVS is the domain and range of the operator.
  Operator_Edge(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                Teuchos::ParameterList& plist) :
      Operator(cvs, plist, OPERATOR_SCHEMA_DOFS_EDGE) {
    set_schema_string("EDGE");
  }

  // rhs update which multiplies by cell
  virtual void UpdateRHS(const CompositeVector& source, bool volume_included);

  // visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Cell_Edge& op,
          const CompositeVector& X, CompositeVector& Y) const;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const;
  
  // visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Cell_Edge& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

    

