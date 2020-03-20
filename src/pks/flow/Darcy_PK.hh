/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1) 
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DARCY_PK_HH_
#define AMANZI_DARCY_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "FnBaseDefs.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Diffusion.hh"
#include "PK_Factory.hh"
#include "TimestepController.hh"
#include "TreeVector.hh"

// Amanzi::Flow
#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

class Darcy_PK : public Flow_PK {
 public:
  Darcy_PK(Teuchos::ParameterList& pk_tree,
           const Teuchos::RCP<Teuchos::ParameterList>& glist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& soln);

  Darcy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
           const std::string& pk_list_name,
           Teuchos::RCP<State> S,
           const Teuchos::RCP<TreeVector>& soln);

  ~Darcy_PK();

  // methods required for PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual void set_dt(double dt) override { dt_ = dt; dt_desirable_ = dt; }
  virtual double get_dt() override { return dt_desirable_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) override; 
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override;

  virtual std::string name() override { return "darcy"; }

  // methods required for time integration interface. It is not used by simple Darcy flow,
  // however, coupled method may need the residual evaluation routine.
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void FunctionalResidual(const double t_old, double t_new, 
                  Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                  Teuchos::RCP<TreeVector> f) override;
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  // -- preconditioner management
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt) override;
 
  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> u) override { return false; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  virtual bool ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override { return false; }

  // -- possibly modifies the correction, after the nonlinear solver (i.e., NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) override {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  // -- experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in state.
  virtual void ChangedSolution() override {};

  // other members of the PK linear solvers
  void SolveFullySaturatedProblem(CompositeVector& u, bool wells_on);

  // access methods
  Teuchos::RCP<Operators::Operator> op() { return op_; }
  Teuchos::RCP<Operators::PDE_Diffusion> op_diff() { return op_diff_; }

 private:
  void InitializeFields_();
  void UpdateSpecificYield_();
  double ErrorEstimate_(double* dTfactor);
  void InitializeStatistics_(bool init_darcy);

  // support of coupled PKs
  void UpdateMatrixBCsUsingFracture_();
  void UpdateSourceUsingMatrix_();
  void FractureConservationLaw_();
  
 protected:
  Teuchos::RCP<TreeVector> soln_;
  
 private:
  Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::PDE_Diffusion> op_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;

  int error_control_;
  double dt_desirable_, dt_factor_;
  std::vector<std::pair<double, double> > dt_history_;  // statistics

  std::string solver_name_;
  bool initialize_with_darcy_;
  int num_itrs_;

  bool flow_on_manifold_;  // true for the DFN model
  bool coupled_to_matrix_, coupled_to_fracture_;

  Teuchos::RCP<CompositeVector> solution;  // next pressure state
  Teuchos::RCP<Epetra_Vector> pdot_cells_prev;  // time derivative of pressure
  Teuchos::RCP<Epetra_Vector> pdot_cells;
  Teuchos::RCP<TimestepController> ts_control_;

  Teuchos::RCP<CompositeVector> specific_yield_copy_;

 private:
  // factory registration
  static RegisteredPKFactory<Darcy_PK> reg_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif

