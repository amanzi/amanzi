/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are EDGEs.
  The only thing really implemented here is the visitor pattern Op
  acceptors. Everything else should be done in the base class, with 
  the exception of special assembly issues.
*/

#ifndef AMANZI_OPERATOR_EDGE_HH_
#define AMANZI_OPERATOR_EDGE_HH_

#include "DenseMatrix.hh"
#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Operator_Edge : public Operator {
 public:
  // The input CVS is the domain and range of the operator.
  Operator_Edge(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                Teuchos::ParameterList& plist) :
      Operator(cvs, plist, OPERATOR_SCHEMA_DOFS_EDGE) {
    set_schema_string("EDGE");
    cell_max_edges = AmanziMesh::MeshAlgorithms::getMaxCellNumEdges(*mesh_);
  }

  // copy constructor
  virtual Teuchos::RCP<Operator> Clone() const;

  // rhs update which multiplies by cell
  virtual void UpdateRHS(const CompositeVector& source, bool volume_included);

  // visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Cell_Edge& op,
          const CompositeVector& X, CompositeVector& Y) const;

  virtual int ApplyMatrixFreeOp(const Op_Edge_Edge& op,
          const CompositeVector& X, CompositeVector& Y) const;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const;
  
  virtual void SymbolicAssembleMatrixOp(const Op_Edge_Edge& op,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const;

  // visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Cell_Edge& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const;

  virtual void AssembleMatrixOp(const Op_Edge_Edge& op,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const;

 protected:
  int cell_max_edges;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

    

