/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)

  Helper class for discrete PDE operators. It provides support of
  common functionality.
*/

#ifndef AMANZI_OPERATOR_PDE_HELPER_HH_
#define AMANZI_OPERATOR_PDE_HELPER_HH_

#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "BCsList.hh"
#include "Mesh.hh"
#include "Op.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_Helper : public BCsList {
 public:
  PDE_Helper() {};
  PDE_Helper(const Teuchos::RCP<Operator>& global_op);
  PDE_Helper(const Teuchos::RCP<AmanziMesh::Mesh>& mesh);
  PDE_Helper(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~PDE_Helper() {};

  // boundary conditions
  virtual void ApplyBCs(bool primary, bool eliminate);

  // access
  // -- global operator (collection of ops with Apply, etc)
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  // -- local operator (container of elemental matrices)
  Teuchos::RCP<Op> local_operator() { return local_op_; }

 protected:
  void ApplyBCs_Cell_Scalar_(const Teuchos::Ptr<BCs>& bc, Teuchos::RCP<Op> op,
                             bool primary, bool eliminate);
  
  void ApplyBCs_Cell_Point_(const Teuchos::Ptr<BCs>& bc, Teuchos::RCP<Op> op,
                            bool primary, bool eliminate);

 private:
  void PopulateDimensions_();

 protected:
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  OperatorType operator_type_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;
  int nedges_owned, nedges_wghost;

  // discretization method
  SpaceNickName space_col_, space_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


