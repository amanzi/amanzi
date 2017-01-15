/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)

  Helper class for discrete PDE operators.
*/

#ifndef AMANZI_OPERATOR_PDE_HELPER_HH_
#define AMANZI_OPERATOR_PDE_HELPER_HH_

#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "Mesh.hh"
#include "Op.hh"
#include "Operator.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_Helper {
 public:
  PDE_Helper() {};
  PDE_Helper(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~PDE_Helper() {};

  void ApplyBCs_Node(const Teuchos::Ptr<BCs>& bc, Teuchos::RCP<Op> op,
                     bool primary, bool eliminate);

  void ApplyBCs_Face(const Teuchos::Ptr<BCs>& bc, Teuchos::RCP<Op> op,
                     bool primary, bool eliminate);
  
 protected:
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  OperatorType operator_type_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


