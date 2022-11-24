/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Special purpose operator, this takes FaceCell Ops and assembles
  a matrix with only the Face entries, enabling the solution of
  consistent faces.
*/

#ifndef AMANZI_OPERATOR_WITH_CONSISTENT_FACE_HH_
#define AMANZI_OPERATOR_WITH_CONSISTENT_FACE_HH_

#include "DenseMatrix.hh"
#include "OperatorDefs.hh"
#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Operator_ConsistentFace : public Operator {
 public:
  // The CVS is the domain and range of the operator
  Operator_ConsistentFace(const Teuchos::RCP<const CompositeVectorSpace>& cvs,
                          Teuchos::ParameterList& plist)
    : Operator(cvs, plist, OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE)
  {
    cell_max_faces_ = mesh_->cell_get_max_faces();
  }

  // copy constructor
  virtual Teuchos::RCP<Operator> Clone() const override;

  // cannot do an assembled forward apply as the assembled thing is not the full
  // operator
  virtual int Apply(const CompositeVector& X, CompositeVector& Y) const override
  {
    return Apply(X, Y, 0.0);
  }
  virtual int Apply(const CompositeVector& X, CompositeVector& Y, double scalar) const override
  {
    return ApplyUnassembled(X, Y, scalar);
  }

  // visit methods for Apply
  virtual int ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                const CompositeVector& X,
                                CompositeVector& Y) const override;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const override;

  // visit methods for assemble
  virtual void AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const override;

 protected:
  int cell_max_faces_;
};

} // namespace Operators
} // namespace Amanzi

#endif
