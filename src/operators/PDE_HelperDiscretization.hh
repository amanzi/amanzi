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
#include "OperatorDefs.hh"
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
  // -- modify matrix due to boundary conditions 
  virtual void ApplyBCs(bool primary, bool eliminate);

  // postprocessing
  // -- flux calculation uses potential p to calculate flux u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) = 0;

  // access
  // -- global operator (collection of ops with Apply, etc)
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  void set_global_operator(const Teuchos::RCP<Operator>& global_op) { global_op_ = global_op; }

  // -- local operator (container of elemental matrices)
  Teuchos::RCP<Op> local_operator() { return local_op_; }
  void set_local_operator(const Teuchos::RCP<Op>& op) { local_op_ = op; }

 protected:
  void ApplyBCs_Cell_Scalar_(const BCs& bc, Teuchos::RCP<Op> op,
                             bool primary, bool eliminate);
  
  void ApplyBCs_Cell_Point_(const BCs& bc, Teuchos::RCP<Op> op,
                            bool primary, bool eliminate);

  void ApplyBCs_Cell_Vector_(const BCs& bc, Teuchos::RCP<Op> op,
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


