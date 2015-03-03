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

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BCs.hh"
#include "BDF1_TI.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"
#include "Upwind.hh"

#include "Flow_PK.hh"
#include "RelativePermeability.hh"
#include "RelPermEvaluator.hh"

#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const std::string& pk_list_name,
              Teuchos::RCP<State> S);
  ~Richards_PK();

  // methods required for PK interface
  void Setup();
  void Initialize();
  bool Advance(double dT_MPC, double& dT_actual); 

  double get_dt() { return dT; }
  void set_dt(double dTnew) { dT = dTnew; dT_desirable_ = dTnew; }

  void CommitState(double dt, const Teuchos::Ptr<State>& S);
  void CalculateDiagnostics(const Teuchos::Ptr<State>& S);

  // main flow methods
  void InitTimeInterval();
  void InitializeAuxiliaryData();

  void UpdateSourceBoundaryData(double T0, double T1, const CompositeVector& u);
  double ErrorNormSTOMP(const CompositeVector& u, const CompositeVector& du);
 
  // methods required for time integration interface
  void Functional(const double T0, double T1, 
                  Teuchos::RCP<CompositeVector> u_old, Teuchos::RCP<CompositeVector> u_new, 
                  Teuchos::RCP<CompositeVector> f);
  void ApplyPreconditioner(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> Hu);
  void UpdatePreconditioner(double T, Teuchos::RCP<const CompositeVector> u, double dT);
  double ErrorNorm(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<const CompositeVector> du);
  bool IsAdmissible(Teuchos::RCP<const CompositeVector> up) { return true; }
  bool ModifyPredictor(double dT, Teuchos::RCP<const CompositeVector> u0,
                       Teuchos::RCP<CompositeVector> u);
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dT, Teuchos::RCP<const CompositeVector> res,
                       Teuchos::RCP<const CompositeVector> u, 
                       Teuchos::RCP<CompositeVector> du);
  void ChangedSolution() {};

  // initization members
  void SolveFullySaturatedProblem(double T0, CompositeVector& u, const std::string& solver_name);
  void EnforceConstraints(double T1, CompositeVector& u);

  void ClipHydrostaticPressure(const double pmin, Epetra_MultiVector& p);
  void ClipHydrostaticPressure(const double pmin, const double s0, Epetra_MultiVector& p);

  int AdvanceToSteadyState_Picard(Teuchos::ParameterList& picard_list);
  double CalculateRelaxationFactor(const Epetra_MultiVector& uold, const Epetra_MultiVector& unew);

  // access methods
  Teuchos::RCP<Operators::Operator> op_matrix() { return op_matrix_; }
  const Teuchos::RCP<CompositeVector> get_solution() { return solution; }
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pressure_eval() { return pressure_eval_; }
  Teuchos::RCP<BDF1_TI<CompositeVector, CompositeVectorSpace> > get_bdf1_dae() { return bdf1_dae; }

  // developement members
  template <class Model> 
  double DeriveBoundaryFaceValue(int f, const CompositeVector& u, const Model& model);
  virtual double BoundaryFaceValue(int f, const CompositeVector& pressure);
  void PlotWRMcurves(Teuchos::ParameterList& plist);

 private:
  void InitializeUpwind_();

  void Functional_AddVaporDiffusion_(Teuchos::RCP<CompositeVector> f);
  void CalculateVaporDiffusionTensor_();

 private:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> rp_list_;

  Teuchos::RCP<RelativePermeability> rel_perm_;

  // solvers
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_, op_pc_solver_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_matrix_diff_, op_preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_acc_;
  Teuchos::RCP<Operators::Upwind<RelativePermeability> > upwind_;
  Teuchos::RCP<Operators::BCs> op_bc_;
  std::string preconditioner_name_, solver_name_, solver_name_constraint_;

  // coupling with energy
  Teuchos::RCP<Operators::Operator> op_vapor_matrix_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_vapor_matrix_diff_;
  bool vapor_diffusion_;
  std::vector<WhetStone::Tensor> K_vapor; 

  // time integrators
  Teuchos::RCP<BDF1_TI<CompositeVector, CompositeVectorSpace> > bdf1_dae;
  int error_control_, num_itrs_;
  double dT_desirable_;
  std::vector<std::pair<double, double> > dT_history_;
  bool initialize_with_darcy_;  // global state of initialization.

  Teuchos::RCP<Epetra_Vector> pdot_cells_prev;  // time derivative of pressure
  Teuchos::RCP<Epetra_Vector> pdot_cells;

  double functional_max_norm;
  int functional_max_cell;

  // water retention models
  Teuchos::RCP<WRMPartition> wrm_;

  // copies of state variables
  Teuchos::RCP<CompositeVector> solution;
  Teuchos::RCP<CompositeVector> darcy_flux_copy;

  // upwind
  int update_upwind;
  Teuchos::RCP<CompositeVector> darcy_flux_upwind;

  // evaluators
  Teuchos::RCP<RelPermEvaluator> rel_perm_eval_;

 private:
  void operator=(const Richards_PK& RPK);

  friend class Richards_PK_Wrapper;
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

