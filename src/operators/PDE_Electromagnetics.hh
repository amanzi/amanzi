/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_ELECTROMAGNETICS_HH_
#define AMANZI_OPERATOR_PDE_ELECTROMAGNETICS_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

class PDE_Electromagnetics : public PDE_HelperDiscretization {
 public:
  PDE_Electromagnetics(const Teuchos::RCP<Operator>& global_op) :
      PDE_HelperDiscretization(global_op),
      K_(Teuchos::null) {};

  PDE_Electromagnetics(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_HelperDiscretization(mesh),
      K_(Teuchos::null)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_ELECTROMAGNETICS;
    Init_(plist);
  }

  virtual ~PDE_Electromagnetics() = default;
  
  // main virtual members
  // -- setup 
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);

  // -- creation of a linearized operator
  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p) override;
  // -- modify matrix due to boundary conditions 
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial 
  //      function, i.e. zeros go in the corresponding matrix columns and 
  //      right-hand side is modified using BC values. This is the optional 
  //      parameter that enforces symmetry for a symmetric tree  operators.
  //    essential_eqn=true indicates that the operator places a positive number on 
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementtion trick.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

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
  int schema_dofs() { return local_op_schema_; }

 protected:
  void Init_(Teuchos::ParameterList& plist);
  void ApplyBCs_Edge_(const Teuchos::Ptr<const BCs>& bc_f,
                      const Teuchos::Ptr<const BCs>& bc_e,
                      bool primary, bool eliminate, bool essential_eqn);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  bool K_symmetric_;

  // operator
  int global_op_schema_, local_op_schema_;

  // miscaleneous
  int mfd_primary_, mfd_secondary_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


