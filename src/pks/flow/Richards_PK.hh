/*
  This is the flow component of the Amanzi code. 

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
#include "BCs.hh"
#include "BDF1_TI.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"
#include "Upwind.hh"

// Flow
#include "Flow_PK.hh"
#include "RelPerm.hh"
#include "RelPermEvaluator.hh"
#include "WRMPartition.hh"

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

  ~Richards_PK();

  // methods required for PK interface
  virtual void Setup();
  virtual void Initialize();

  virtual double get_dt() { return dt_; }
  virtual void set_dt(double dt) { dt_ = dt; dt_desirable_ = dt_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics();

  virtual std::string name() { return passwd_; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void Functional(const double t_old, double t_new, 
                  Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                  Teuchos::RCP<TreeVector> f);
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // -- management of the preconditioner
  void ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu);
  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt);

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) { return true; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  bool ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u);

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, 
                       Teuchos::RCP<TreeVector> du);

  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  void ChangedSolution() {
    pressure_eval_->SetFieldAsChanged(S_.ptr());
  }

  // other flow methods
  // -- initization members
  void SolveFullySaturatedProblem(double t_old, CompositeVector& u, const std::string& solver_name);
  void EnforceConstraints(double t_new, Teuchos::RCP<CompositeVector> u);
  void EnforceConstraints_Inflow(Teuchos::RCP<CompositeVector> u);

  void ClipHydrostaticPressure(const double pmin, Epetra_MultiVector& p);
  void ClipHydrostaticPressure(const double pmin, const double s0, Epetra_MultiVector& p);

  int AdvanceToSteadyState_Picard(Teuchos::ParameterList& picard_list);
  double CalculateRelaxationFactor(const Epetra_MultiVector& uold, const Epetra_MultiVector& unew);

  // -- miscaleneous methods
  void UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u);
  double ErrorNormSTOMP(const CompositeVector& u, const CompositeVector& du);
 
  // access methods
  Teuchos::RCP<Operators::Operator> op_matrix() { return op_matrix_; }
  const Teuchos::RCP<CompositeVector> get_solution() { return solution; }
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > get_bdf1_dae() { return bdf1_dae; }

  // developement members
  template <class Model> 
  double DeriveBoundaryFaceValue(int f, const CompositeVector& u, const Model& model);
  virtual double BoundaryFaceValue(int f, const CompositeVector& pressure);
  void PlotWRMcurves(Teuchos::ParameterList& plist);

 private:
  void InitializeFields_();
  void InitializeUpwind_();

  void Functional_AddVaporDiffusion_(Teuchos::RCP<CompositeVector> f);
  void CalculateVaporDiffusionTensor_(Teuchos::RCP<CompositeVector>& kvapor_pres,
                                      Teuchos::RCP<CompositeVector>& kvapor_temp);

 private:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> rp_list_;

  // pointerds to primary field
  const Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution;

  // water retention models
  Teuchos::RCP<WRMPartition> wrm_;

  Teuchos::RCP<RelPerm> relperm_;
  int krel_upwind_method_;
  Teuchos::RCP<CompositeVector> krel_;
  Teuchos::RCP<CompositeVector> dKdP_;

  // solvers
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_, op_pc_solver_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_matrix_diff_, op_preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_acc_;
  Teuchos::RCP<Operators::Upwind<RelPerm> > upwind_;
  Teuchos::RCP<Operators::BCs> op_bc_;
  std::string preconditioner_name_, solver_name_, solver_name_constraint_;

  // coupling with energy
  Teuchos::RCP<Operators::Operator> op_vapor_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_vapor_diff_;
  bool vapor_diffusion_;

  // time integrators
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae;
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
  int update_upwind;
  Teuchos::RCP<CompositeVector> darcy_flux_upwind;

  // evaluators
  Teuchos::RCP<RelPermEvaluator> rel_perm_eval_;

 private:
  void operator=(const Richards_PK& RPK);

 private:
  // factory registration
  static RegisteredPKFactory<Richards_PK> reg_;
};


/* ******************************************************************
* Calculates solution value on the boundary.
****************************************************************** */
template <class Model> 
double Richards_PK::DeriveBoundaryFaceValue(
    int f, const CompositeVector& u, const Model& model) 
{
  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    return u_face[f][0];
  } else {
    const std::vector<int>& bc_model = op_bc_->bc_model();
    const std::vector<double>& bc_value = op_bc_->bc_value();

    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      return bc_value[f];
    } else {
      const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];
      return u_cell[0][c];
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

#endif

