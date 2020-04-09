/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

//! Composes BCs, local_op, and global_operator objects for use by discretizations.


#ifndef AMANZI_OPERATOR_PDE_HELPER_DISCRETIZATION_HH_
#define AMANZI_OPERATOR_PDE_HELPER_DISCRETIZATION_HH_

#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "Mesh.hh"
#include "Op.hh"
#include "Operator.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_HelperDiscretization {
 public:
  PDE_HelperDiscretization(const Teuchos::RCP<Operator>& global_op);
  PDE_HelperDiscretization(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // access
  // -- global operator (collection of ops with Apply, etc)
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  void set_global_operator(const Teuchos::RCP<Operator>& global_op) { global_op_ = global_op; }

  // -- local operator (container of elemental matrices)
  Teuchos::RCP<Op> local_op() { return local_op_; }
  Teuchos::RCP<const Op> local_op() const { return local_op_; }
  void set_local_op(const Teuchos::RCP<Op>& op);

  // boundary conditions (BC) require information on test and
  // trial spaces. For a single PDE, these BCs could be the same.
  void SetBCs(const Teuchos::RCP<const BCs>& bc_trial,
              const Teuchos::RCP<const BCs>& bc_test);
  void AddBCs(const Teuchos::RCP<const BCs>& bc_trial,
              const Teuchos::RCP<const BCs>& bc_test);

 private:
  void PopulateDimensions_();

 protected:
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  std::vector<Teuchos::RCP<const BCs> > bcs_trial_, bcs_test_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;
  int nedges_owned, nedges_wghost;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


