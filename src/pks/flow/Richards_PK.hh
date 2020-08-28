/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_RICHARDS_PK_HH_
#define AMANZI_RICHARDS_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BDF1_TI.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Diffusion.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"
#include "Upwind.hh"

// Amanzi::Flow
#include "Flow_PK.hh"
#include "MultiscaleFlowPorosityPartition.hh"
#include "RelPerm.hh"
#include "RelPermEvaluator.hh"
#include "WRMPartition.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln);

  Richards_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const std::string& pk_list_name,
              Teuchos::RCP<State> S,
              const Teuchos::RCP<TreeVector>& soln);

  ~Richards_PK() {};

  // methods required for PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual double get_dt() override { return dt_; }
  virtual void set_dt(double dt) override { dt_ = dt; dt_desirable_ = dt_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) override;
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override;

  virtual std::string name() override { return "richards"; }

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
                               Teuchos::RCP<TreeVector> u) override;

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, 
                       Teuchos::RCP<TreeVector> du) override;

  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  virtual void ChangedSolution() override {
    pressure_eval_->SetFieldAsChanged(S_.ptr());
  }

  // -- returns the number of linear iterations.
  virtual int ReportStatistics() override { 
    return op_preconditioner_->apply_calls();
  }

  // other flow methods
  // -- initialization members
  void SolveFullySaturatedProblem(double t_old, CompositeVector& u, const std::string& solver_name);
  void EnforceConstraints(double t_new, Teuchos::RCP<CompositeVector> u);
  void UpwindInflowBoundary(Teuchos::RCP<const CompositeVector> u);
  void UpwindInflowBoundary_New(Teuchos::RCP<const CompositeVector> u);

  void ClipHydrostaticPressure(const double pmin, Epetra_MultiVector& p);
  void ClipHydrostaticPressure(const double pmin, const double s0, Epetra_MultiVector& p);

  int AdvanceToSteadyState_Picard(Teuchos::ParameterList& picard_list);
  double CalculateRelaxationFactor(const Epetra_MultiVector& uold, const Epetra_MultiVector& unew);

  // -- mutiscale methods
  void CalculateVWContentMatrix_();
  void VV_ReportMultiscale();

  // -- miscaleneous methods
  double ErrorNormSTOMP(const CompositeVector& u, const CompositeVector& du);

  // -- access methods
  virtual Teuchos::RCP<Operators::Operator>
      my_operator(const Operators::OperatorType& type) override;

  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization>
      my_pde(const Operators::PDEType& type) override { return op_matrix_diff_; } 

  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > get_bdf1_dae() { return bdf1_dae_; }

  // -- verbose output and visualization methods
  void PlotWRMcurves(Teuchos::ParameterList& plist);

  // -- developement methods
  double DeriveBoundaryFaceValue(int f, const CompositeVector& u, Teuchos::RCP<const WRM> model);
  virtual double BoundaryFaceValue(int f, const CompositeVector& pressure) override;

 private:
  void InitializeFields_();
  void InitializeFieldFromField_(const std::string& field0, const std::string& field1, bool call_evaluator);
  void InitializeStatistics_();

  void Functional_AddVaporDiffusion_(Teuchos::RCP<CompositeVector> f);
  void CalculateVaporDiffusionTensor_(Teuchos::RCP<CompositeVector>& kvapor_pres,
                                      Teuchos::RCP<CompositeVector>& kvapor_temp);

  void Functional_AddMassTransferMatrix_(double dt, Teuchos::RCP<CompositeVector> f);

  // The water content change in a cell equals exactly to the balance of Darcy fluxes.
  // This balance leads to a monotone translport.
  void CalculateCNLSLimiter_(const CompositeVector& wc, const CompositeVector& dwc_dp, double tol);
  void ApplyCNLSLimiter_();
 
 private:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;

  // pointerds to primary field
  const Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution;

  // water retention models
  Teuchos::RCP<WRMPartition> wrm_;

  Teuchos::RCP<RelPerm> relperm_;
  Teuchos::RCP<CompositeVector> krel_;
  Teuchos::RCP<CompositeVector> dKdP_;

  // solvers
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_;
  Teuchos::RCP<Matrix<CompositeVector,CompositeVectorSpace>> op_pc_solver_;
  Teuchos::RCP<Operators::PDE_Diffusion> op_matrix_diff_, op_preconditioner_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;
  Teuchos::RCP<Operators::Upwind<RelPerm> > upwind_;
  std::string solver_name_, solver_name_constraint_;

  // coupling with energy
  Teuchos::RCP<Operators::Operator> op_vapor_;
  Teuchos::RCP<Operators::PDE_Diffusion> op_vapor_diff_;
  bool vapor_diffusion_;

  // multiscale models
  bool multiscale_porosity_;
  int ms_itrs_, ms_calls_;
  Teuchos::RCP<MultiscaleFlowPorosityPartition> msp_;

  // time integrators
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae_;
  int error_control_, num_itrs_;
  double dt_desirable_;
  std::vector<std::pair<double, double> > dT_history_;
  bool initialize_with_darcy_;  // global state of initialization.

  Teuchos::RCP<Epetra_Vector> pdot_cells_prev;  // time derivative of pressure
  Teuchos::RCP<Epetra_Vector> pdot_cells;

  double functional_max_norm;
  int functional_max_cell;

  // copies of state fields
  Teuchos::RCP<CompositeVector> darcy_flux_copy;

  // upwind
  int upwind_frequency_;

  // evaluators
  Teuchos::RCP<RelPermEvaluator> rel_perm_eval_;

  // consistent water content and Darcy fluxes
  bool algebraic_water_content_balance_;
  Teuchos::RCP<CompositeVector> cnls_limiter_;

 private:
  void operator=(const Richards_PK& RPK);

 private:
  // factory registration
  static RegisteredPKFactory<Richards_PK> reg_;
};



}  // namespace Flow
}  // namespace Amanzi

#endif

