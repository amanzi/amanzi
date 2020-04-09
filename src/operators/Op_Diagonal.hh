/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
          Konstantin Lipnikov (lipnikov@lanl.gov)

  This operator is a container for local matrices of variable
  length placed on operator's diagonal. The locations are
  specified by two sets of indices.
*/

#ifndef AMANZI_OP_DIAGONAL_HH_
#define AMANZI_OP_DIAGONAL_HH_

#include <memory>
#include <string>
#include <vector>

#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Op_Diagonal : public Op {
 public:
  Op_Diagonal(const std::string& name,
              std::string row_compname, std::string col_compname,
              std::shared_ptr<const std::vector<std::vector<int> > >& row_inds,
              std::shared_ptr<const std::vector<std::vector<int> > >& col_inds) :
      Op(OPERATOR_SCHEMA_INDICES, name),
      row_compname_(row_compname),
      col_compname_(col_compname),
      row_inds_(row_inds),
      col_inds_(col_inds) {
    WhetStone::DenseMatrix null_matrix;
    matrices.resize(row_inds->size(), null_matrix);
    matrices_shadow = matrices;
  }
  ~Op_Diagonal() {};

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const {
    assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }


  // incomplete members
  virtual void Rescale(const CompositeVector& scaling) { AMANZI_ASSERT(0); } 

  // access 
  const std::vector<std::vector<int> >& row_inds() const { return *row_inds_; }
  const std::vector<std::vector<int> >& col_inds() const { return *col_inds_; }
  
  std::string row_compname() const { return row_compname_; }
  std::string col_compname() const { return col_compname_; }

 private:
  std::string row_compname_, col_compname_;
  std::shared_ptr<const std::vector<std::vector<int> > > row_inds_, col_inds_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


