/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Examples of usage this operator are in test/operators_elasticity.cc
  and test/operators_stokes.cc
*/

#ifndef AMANZI_OPERATOR_PDE_ELASTICITY_HH_
#define AMANZI_OPERATOR_PDE_ELASTICITY_HH_

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

// Amanzi::Operators
#include "BilinearForm.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class PDE_Elasticity : public PDE_HelperDiscretization {
 public:
  PDE_Elasticity(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_HelperDiscretization(mesh),
      K_(Teuchos::null),
      K_default_(1.0)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_ELASTICITY;
    Init_(plist);
  }

  // main virtual members
  // -- setup 
  void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
  void SetTensorCoefficient(double K);

  // -- creation of an operator
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;

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
  //      implementtion trick/
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

  // -- postprocessing: calculated stress u from displacement p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

  // access
  const Schema& global_schema_col() { return global_schema_col_; }
  const Schema& schema_col() { return local_schema_col_; }
  const Schema& schema_row() { return local_schema_row_; }

 protected:
  void Init_(Teuchos::ParameterList& plist);

 private:
  // this should be generalized when the next use case appear
  void ApplyBCs_Node_Point_(const BCs& bc, Teuchos::RCP<Op> op,
                            bool primary, bool eliminate, bool essential_eqn);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  double K_default_;

  Teuchos::RCP<WhetStone::BilinearForm> mfd_;
  AmanziMesh::Entity_kind base_;

  // operator and schemas
  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


