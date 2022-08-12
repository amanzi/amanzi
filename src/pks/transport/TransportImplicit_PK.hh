/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Implementation of implicit time integration algorithms.
*/

#ifndef AMANZI_TRANSPORT_IMPLICIT_PK_HH_
#define AMANZI_TRANSPORT_IMPLICIT_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseVector.hh"
#include "Key.hh"
#include "LimiterCell.hh"
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"

// Amanzi
#include "Transport_PK.hh"
#include "TransportDefs.hh"
#include "TransportDomainFunction.hh"

#include "BCs.hh"
#include "BDF1_TI.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PK_BDF.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace Transport {

class TransportImplicit_PK : public Transport_PK,
                             public PK_BDF {
 public:
  TransportImplicit_PK(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& glist,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& soln);

  TransportImplicit_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                       Teuchos::RCP<State> S, 
                       const std::string& pk_list_name,
                       std::vector<std::string>& component_names);

  ~TransportImplicit_PK() {};

  virtual void Initialize() override;  
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) override ; 

  virtual double get_dt() override { return dt_; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  virtual void FunctionalResidual(
      const double t_old, double t_new, 
      Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
      Teuchos::RCP<TreeVector> f) override;
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  // -- management of the preconditioner
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt) override;

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override { return true; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  virtual bool ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override { return false; }

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, 
                       Teuchos::RCP<TreeVector> du) override {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  virtual void ChangedSolution() override {};

  void UpdateLinearSystem(double t_old, double t_new, int component);
  void UpdateBoundaryData(double t_old, double t_new, int component);

  // access
  Teuchos::RCP<Operators::Operator> op() { return op_; }
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_adv() { return op_adv_; }
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc() { return op_acc_; }

 private:  
  bool AdvanceStepLO_(double t_old, double t_new, int* tot_itrs);
  bool AdvanceStepHO_(double t_old, double t_new, int* tot_itrs);

 private:
  // solvers
  Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_adv_;
  Teuchos::RCP<Operators::PDE_Diffusion> op_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;
  Teuchos::RCP<Operators::BCs> op_bc_;
  std::string solver_name_, solver_name_constraint_;

  Teuchos::RCP<CompositeVector> solution_;
  std::vector<Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>>> bdf1_dae_;
  Teuchos::RCP<Matrix<CompositeVector,CompositeVectorSpace>> op_pc_solver_;

  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  // factory registration
  static RegisteredPKFactory<TransportImplicit_PK> reg_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif
