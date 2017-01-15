/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_ELECTROMAGNETICS_HH_
#define AMANZI_OPERATOR_ELECTROMAGNETICS_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

#include "BCsList.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class Electromagnetics : public BCsList {
 public:
  Electromagnetics(const Teuchos::RCP<Operator>& global_op) :
      global_op_(global_op),
      K_(Teuchos::null),
      ncells_owned(-1),
      ncells_wghost(-1),
      nfaces_owned(-1),
      nfaces_wghost(-1),
      nedges_owned(-1),
      nedges_wghost(-1)
  {};

  Electromagnetics(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      plist_(plist),
      mesh_(mesh),
      K_(Teuchos::null),
      global_op_(Teuchos::null),
      ncells_owned(-1),
      ncells_wghost(-1),
      nfaces_owned(-1),
      nfaces_wghost(-1),
      nedges_owned(-1),
      nedges_wghost(-1)
  {
    operator_type_ = OPERATOR_ELECTROMAGNETICS;
    InitElectromagnetics_(plist);
  }

  // main virtual members
  // -- setup 
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);

  // -- creation of an operator
  virtual void UpdateMatrices();

  // -- matrix modification
  virtual void ApplyBCs(bool primary, bool eliminate);

  // new virtual members
  // -- before solving the problem
  virtual void ModifyMatrices(CompositeVector& E, CompositeVector& B, double dt) {};

  // -- after solving the problem
  virtual void ModifyFields(CompositeVector& E, CompositeVector& B, double dt) {};

  // access
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  int schema_prec_dofs() { return global_op_schema_; }

  Teuchos::RCP<const Op> local_matrices() const { return local_op_; }
  Teuchos::RCP<Op> local_matrices() { return local_op_; }
  int schema_dofs() { return local_op_schema_; }

 protected:
  void InitElectromagnetics_(Teuchos::ParameterList& plist);
  void ApplyBCs_Edge_(const Teuchos::Ptr<BCs>& bc_f,
                      const Teuchos::Ptr<BCs>& bc_e, bool primary, bool eliminate);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  bool K_symmetric_;

  // operator
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  int global_op_schema_, local_op_schema_;
  OperatorType operator_type_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nedges_owned, nedges_wghost;

  // miscaleneous
  Teuchos::ParameterList plist_;
  int mfd_primary_, mfd_secondary_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


