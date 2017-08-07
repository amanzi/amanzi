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

#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Helper.hh"

namespace Amanzi {
namespace Operators {

class Electromagnetics : public PDE_Helper {
 public:
  Electromagnetics(const Teuchos::RCP<Operator>& global_op) :
      PDE_Helper(global_op),
      K_(Teuchos::null) {};

  Electromagnetics(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Helper(mesh),
      plist_(plist),
      K_(Teuchos::null)
  {
    global_op_ = Teuchos::null;
    operator_type_ = OPERATOR_ELECTROMAGNETICS;
    Init_(plist);
  }

  virtual ~Electromagnetics() = default;
  
  // main virtual members
  // -- setup 
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
  // -- creation of a linearized operator
  using PDE_Helper::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;
  // -- modification of the operator due to boundary conditions
  //    Variable 'primary' indicates that we put 1 on the matrix diagonal.
  //    Variable 'eliminate' says that we eliminate essential BCs for the 
  //    trial function, i.e. zeros go in the corresponding matrix columns.
  virtual void ApplyBCs(bool primary, bool eliminate);

  // -- postprocessing: calculated flux u from potential p
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

  // new virtual members
  // -- before solving the problem
  virtual void ModifyMatrices(CompositeVector& E, CompositeVector& B, double dt) {};

  // -- after solving the problem
  virtual void ModifyFields(CompositeVector& E, CompositeVector& B, double dt) {};

  // access
  int schema_prec_dofs() { return global_op_schema_; }

  Teuchos::RCP<const Op> local_matrices() const { return local_op_; }
  Teuchos::RCP<Op> local_matrices() { return local_op_; }
  int schema_dofs() { return local_op_schema_; }

 protected:
  void Init_(Teuchos::ParameterList& plist);
  void ApplyBCs_Edge_(const Teuchos::Ptr<BCs>& bc_f,
                      const Teuchos::Ptr<BCs>& bc_e, bool primary, bool eliminate);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  bool K_symmetric_;

  // operator
  int global_op_schema_, local_op_schema_;

  // miscaleneous
  Teuchos::ParameterList plist_;
  int mfd_primary_, mfd_secondary_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


