/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Operator whose unknowns are CELL + FACE but possibly from multiple
  meshes in contact. Serial implementation.
*/

#ifndef AMANZI_OPERATOR_MULTI_MESH_HH_
#define AMANZI_OPERATOR_MULTI_MESH_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"

#include "DenseMatrix.hh"
#include "Op.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

typedef std::map<int, std::map<int, double>> InterfaceData;

class Operator_MultiMesh : public Operator {
 public:
  Operator_MultiMesh(Teuchos::ParameterList& plist,
                     Teuchos::RCP<Operator> global_op,
                     Teuchos::RCP<Op> local_op,
                     std::vector<int>& interface_block,
                     InterfaceData& interface_data);

  Teuchos::RCP<Operator> Clone() const;

  // visit methods for symbolic assemble
  virtual void SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const;

  // visit methods for actual assemble
  virtual void AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const;

 private:
  int max_size_;
  const std::vector<int>& interface_block_;
  InterfaceData& interface_data_;
};

} // namespace Operators
} // namespace Amanzi

#endif
