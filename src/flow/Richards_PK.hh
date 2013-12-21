/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
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
#include "AztecOO.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BDF1_TI.hh"

#include "Flow_PK.hh"
#include "Matrix.hh"
#include "WaterRetentionModel.hh"
#include "RelativePermeability.hh"
#include "TI_Specs.hh"

namespace Amanzi {
namespace AmanziFlow {

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::ParameterList& global_list, Teuchos::RCP<State> S);
  ~Richards_PK();

  // main methods
  void InitPK();
  void InitSteadyState(double T0, double dT0);
  void InitTransient(double T0, double dT0);
  void InitPicard(double T0);
  void InitNextTI(double T0, double dT0, TI_Specs& ti_specs);

  double CalculateFlowDt();
  int Advance(double dT_MPC); 
  int AdvanceToSteadyState(double T0, double dT0);
  void InitializeAuxiliaryData();
  void InitializeSteadySaturated();

  int AdvanceToSteadyState_Picard(TI_Specs& ti_specs);
  int AdvanceToSteadyState_BackwardEuler(TI_Specs& ti_specs);
  int AdvanceToSteadyState_BDF1(TI_Specs& ti_specs);

  void CommitState(Teuchos::RCP<State> S);

  // methods for experimental time integration
  double ErrorNormSTOMP(const CompositeVector& u, const CompositeVector& du);
  double ErrorNormPicardExperimental(const CompositeVector& uold, const CompositeVector& unew);
 
  // methods required for time integration
  void Functional(const double Told, double Tnew, 
                  Teuchos::RCP<CompositeVector> u_old, Teuchos::RCP<CompositeVector> u_new, 
                  Teuchos::RCP<CompositeVector> f);
  void ApplyPreconditioner(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> Hu);
  void UpdatePreconditioner(double T, Teuchos::RCP<const CompositeVector> u, double dT);
  double ErrorNorm(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<const CompositeVector> du);
  void update_norm(double rtol, double atol) {};
  bool is_admissible(Teuchos::RCP<const CompositeVector> up) { 
   return false; 
  }
  bool ModifyPredictor(double h, Teuchos::RCP<CompositeVector> up) {
    return false;
  }
  bool ModifyCorrection(double h, Teuchos::RCP<const CompositeVector> res,
                         Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> du);
  void changed_solution() {};

  // other main methods
  void AddTimeDerivative_MFD(Epetra_Vector& p, double dTp, Matrix_MFD* matrix_operator);
  void AddTimeDerivative_MFDfake(Epetra_Vector& p, double dTp, Matrix_MFD* matrix_operator);

  double ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot);
  void AssembleMatrixMFD(const CompositeVector &u, double Tp);
  void AssemblePreconditionerMFD(const CompositeVector &u, double Tp, double dTp);

  void UpdateSourceBoundaryData(double Tp, const CompositeVector& pressure);
  double AdaptiveTimeStepEstimate(double* dTfactor);

  // linear problems and solvers
  void AssembleSteadyStateMatrix(FlowMatrix* matrix);
  void AssembleSteadyStatePreconditioner(FlowMatrix* preconditioner);
  void SolveFullySaturatedProblem(double T, CompositeVector& u, LinearSolver_Specs& ls_specs);
  void EnforceConstraints(double Tp, CompositeVector& u);

  // water retention models
  void DeriveSaturationFromPressure(const Epetra_MultiVector& p, Epetra_MultiVector& s);
  void DerivePressureFromSaturation(const Epetra_MultiVector& s, Epetra_MultiVector& p);

  // initization members
  void ClipHydrostaticPressure(const double pmin, Epetra_MultiVector& p);
  void ClipHydrostaticPressure(const double pmin, const double s0, Epetra_MultiVector& p);

  double CalculateRelaxationFactor(const Epetra_MultiVector& uold,
                                   const Epetra_MultiVector& unew);

  // control method
  void ResetErrorControl(int error) { error_control_ = error; }
  void ResetParameterList(const Teuchos::ParameterList& rp_list_new) { rp_list_ = rp_list_new; }
  
  // access methods
  Teuchos::RCP<FlowMatrix> matrix() { return matrix_; }
  Teuchos::RCP<FlowMatrix> preconditioner() { return preconditioner_; }
  const Teuchos::RCP<CompositeVector> get_solution() { return solution; }

  // developement members
  bool SetSymmetryProperty();
  void ImproveAlgebraicConsistency(const Epetra_Vector& ws_prev, Epetra_Vector& ws);
  
  // auxilliary data management
  void UpdateAuxilliaryData();

 public:
  Teuchos::ParameterList rp_list_;

 private:
  Teuchos::RCP<FlowMatrix> matrix_;
  Teuchos::RCP<FlowMatrix> preconditioner_;

  Teuchos::RCP<BDF1_TI<CompositeVector, CompositeVectorSpace> > bdf1_dae;  // Time integrators
  int block_picard;

  int error_control_;
  double dT_desirable_;

  double functional_max_norm;
  int functional_max_cell;

  Teuchos::RCP<CompositeVector> solution;
  Teuchos::RCP<Epetra_Vector> pdot_cells_prev;  // time derivative of pressure
  Teuchos::RCP<Epetra_Vector> pdot_cells;

  Teuchos::RCP<RelativePermeability> rel_perm;

  Teuchos::RCP<Epetra_Vector> Transmis_faces;
  Teuchos::RCP<Epetra_Vector> Grav_term_faces;

  bool is_matrix_symmetric;
  double mass_bc, mass_amanzi;

 private:
  void operator=(const Richards_PK& RPK);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

