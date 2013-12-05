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

#include "BDF2_TI.hh"
#include "BDF1_TI.hh"

#include "Flow_PK.hh"
#include "Matrix_MFD.hh"
#include "WaterRetentionModel.hh"
#include "RelativePermeability.hh"
#include "TI_Specs.hh"

// #include "TimerManager.hh"


namespace Amanzi {
namespace AmanziFlow {

class Flow_State;  // forward declarations

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
  int AdvanceToSteadyState_BDF2(TI_Specs& ti_specs);

  void CommitState(Teuchos::RCP<Flow_State> FS);

  // methods for experimental time integration
  double ErrorNormRC1(const Epetra_Vector& u, const Epetra_Vector& du);
  double ErrorNormSTOMP(const Epetra_Vector& u, const Epetra_Vector& du);
  double ErrorNormPicardExperimental(const Epetra_Vector& uold, const Epetra_Vector& unew);
 
  // methods required for time integration
  void fun(double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs, double dT = 0.0);
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_norm(double rtol, double atol) {};
  void update_precon(double T, const Epetra_Vector& u, double dT, int& ierr);
  bool modify_update_step(double h, Epetra_Vector&u, Epetra_Vector& du );
  bool IsPureNewton() const;

  // other main methods
  void AddTimeDerivative_MFD(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix_operator);
  void AddTimeDerivative_MFDfake(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix_operator);

  double ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot);
  void AssembleMatrixMFD(const Epetra_Vector &u, double Tp);
  void AssemblePreconditionerMFD(const Epetra_Vector &u, double Tp, double dTp);
  void ComputeTransmissibilities(Epetra_Vector& Trans_faces, Epetra_Vector& grav_faces);

  void UpdateSourceBoundaryData(double Tp, CompositeVector& pressure);
  double AdaptiveTimeStepEstimate(double* dTfactor);

  // linear problems and solvers
  void AssembleSteadyStateMatrix_MFD(Matrix_MFD* matrix);
  void AssembleSteadyStatePreconditioner_MFD(Matrix_MFD* preconditioner);
  void SolveFullySaturatedProblem(double T, CompositeVector& u, LinearSolver_Specs& ls_specs);
  void EnforceConstraints_MFD(double Tp, Epetra_Vector& u);

  // water retention models
  void DeriveSaturationFromPressure(const Epetra_MultiVector& p, Epetra_MultiVector& s);
  void DerivePressureFromSaturation(const Epetra_MultiVector& s, Epetra_MultiVector& p);

  // initization members
  void ClipHydrostaticPressure(const double pmin, Epetra_Vector& pressure_cells);
  void ClipHydrostaticPressure(const double pmin, const double s0, Epetra_Vector& pressure_cells);

  double CalculateRelaxationFactor(const Epetra_Vector& uold, const Epetra_Vector& unew);

  // control method
  void ResetErrorControl(int error) { error_control_ = error; }
  void ResetParameterList(const Teuchos::ParameterList& rp_list_new) { rp_list_ = rp_list_new; }
  
  // access methods
  Teuchos::RCP<Matrix_MFD> matrix() { return matrix_; }
  Teuchos::RCP<Matrix_MFD> preconditioner() { return preconditioner_; }

  // developement members
  bool SetSymmetryProperty();
  void ImproveAlgebraicConsistency(const Epetra_Vector& flux, 
                                   const Epetra_Vector& ws_prev, Epetra_Vector& ws);
  
  int ApllyPrecInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

  // auxilliary data management
  void UpdateAuxilliaryData();

 public:
  Teuchos::ParameterList rp_list_;

 private:
  Teuchos::RCP<Matrix_MFD> matrix_;
  Teuchos::RCP<Matrix_MFD> preconditioner_;

  BDF2::Dae* bdf2_dae;  // Time integrators
  BDF1Dae* bdf1_dae;
  int block_picard;

  int error_control_;
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

