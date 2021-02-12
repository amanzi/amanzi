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

#ifndef AMANZI_OPERATOR_PDE_HELPER_DISCRETIZATION_HH_
#define AMANZI_OPERATOR_PDE_HELPER_DISCRETIZATION_HH_

#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "Mesh.hh"
#include "Op.hh"
#include "Operator.hh"
#include "PDE_HelperBCsList.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_HelperDiscretization : public PDE_HelperBCsList {
 public:
  PDE_HelperDiscretization() {};
  PDE_HelperDiscretization(const Teuchos::RCP<Operator>& global_op);
  PDE_HelperDiscretization(const Teuchos::RCP<AmanziMesh::Mesh>& mesh);
  PDE_HelperDiscretization(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~PDE_HelperDiscretization() {};

  // generate linearized operator
  // -- generate matrix. We can use parameter to define coefficeints
  //    or/and perform on-a-fly linearization. 
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) = 0;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u) {
    UpdateMatrices(u, Teuchos::null);
  }
  virtual void UpdateMatrices() {
    UpdateMatrices(Teuchos::null, Teuchos::null);
  }
  // -- modify matrix due to boundary conditions: generic implementation 
  //    for PDE classes based on new schema: 
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial 
  //      function, i.e. zeros go in the corresponding matrix columns and 
  //      right-hand side is modified using BC values. This is the optional 
  //      parameter that enforces symmetry for a symmetric tree  operators.
  //    essential_eqn=true indicates that the operator places a positive number on 
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementtion trick.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn);

  // postprocessing
  // -- flux calculation uses potential p to calculate flux u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) = 0;

  // preprocessing
  // -- set essential boundary conditions
  void EnforceBCs(CompositeVector& field);

  // access
  // -- global operator (collection of ops with Apply, etc)
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  void set_global_operator(const Teuchos::RCP<Operator>& global_op) { global_op_ = global_op; }

  // -- local operator (container of elemental matrices)
  Teuchos::RCP<Op> local_op() { return local_op_; }
  Teuchos::RCP<const Op> local_op() const { return local_op_; }
  void set_local_op(const Teuchos::RCP<Op>& op);
  
  // -- schemas
  const Schema& global_schema_col() { return global_schema_col_; }
  const Schema& schema_col() { return local_schema_col_; }
  const Schema& schema_row() { return local_schema_row_; }

 protected:
  void ApplyBCs_Cell_Scalar_(const BCs& bc, Teuchos::RCP<Op> op,
                             bool primary, bool eliminate, bool essential_eqn);
  
  void ApplyBCs_Cell_Point_(const BCs& bc, Teuchos::RCP<Op> op,
                            bool primary, bool eliminate, bool essential_eqn);

  void ApplyBCs_Cell_Vector_(const BCs& bc, Teuchos::RCP<Op> op,
                             bool primary, bool eliminate, bool essential_eqn);

 private:
  void PopulateDimensions_();

 protected:
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  PDEType pde_type_;

  // schemas
  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;
  int nedges_owned, nedges_wghost;
};


// non-member functions
Teuchos::RCP<CompositeVectorSpace> CreateFracturedMatrixCVS(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const Teuchos::RCP<const AmanziMesh::Mesh>& fracture);

Teuchos::RCP<CompositeVectorSpace> CreateNonManifoldCVS(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

}  // namespace Operators
}  // namespace Amanzi

#endif


