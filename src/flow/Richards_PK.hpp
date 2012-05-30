/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __RICHARDS_PK_HPP__
#define __RICHARDS_PK_HPP__

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "AztecOO.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BDF2_Dae.hpp"
#include "BDF1_Dae.hh"

#include "Flow_PK.hpp"
#include "Flow_BC_Factory.hpp"
#include "Matrix_MFD.hpp"
#include "WaterRetentionModel.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Flow_State;  // forward declarations

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::ParameterList& global_list, Teuchos::RCP<Flow_State> FS_MPC);
  ~Richards_PK();

  // main methods
  void InitPK(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL);
  void InitSteadyState(double T0, double dT0);
  void InitTransient(double T0, double dT0);

  double CalculateFlowDt();
  int Advance(double dT_MPC); 
  int AdvanceToSteadyState();

  void CommitState(Teuchos::RCP<Flow_State> FS);
  void CommitStateForTransport(Teuchos::RCP<Flow_State> FS);
  void DeriveDarcyVelocity(const Epetra_Vector& flux, Epetra_MultiVector& velocity);

  int AdvanceSteadyState_Picard();
  int AdvanceSteadyState_BackwardEuler();
  int AdvanceSteadyState_BDF1();
  int AdvanceSteadyState_BDF2();

  // methods for experimental time integration
  int PicardStep(double T, double dT, double& dTnext);
  double ErrorNormRC1(const Epetra_Vector& u, const Epetra_Vector& du);
  double ErrorNormSTOMP(const Epetra_Vector& u, const Epetra_Vector& du);
  double ErrorNormPicardExperimental(const Epetra_Vector& uold, const Epetra_Vector& unew);
 
  // methods required for time integration
  void fun(double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs, double dT = 0.0);
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_precon(double T, const Epetra_Vector& u, double dT, int& ierr);
  void update_norm(double rtol, double atol) {};

  // other main methods
  void SetAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K);
  void CalculateKVectorUnit(const AmanziGeometry::Point& g, std::vector<AmanziGeometry::Point>& Kg_unit);

  void CalculateRelativePermeability(const Epetra_Vector& u);
  void CalculateRelativePermeabilityCell(const Epetra_Vector& p);
  void CalculateRelativePermeabilityFace(const Epetra_Vector& p);
  void CalculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p);
  void CalculateRelativePermeabilityUpwindFlux(const Epetra_Vector& p, const Epetra_Vector& flux);
  void CalculateRelativePermeabilityArithmeticMean(const Epetra_Vector& p);

  void AddTimeDerivative_MFD(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix);
  void AddTimeDerivative_MFDfake(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix);
  void AddTimeDerivative_MFDpicard(Epetra_Vector& pressure_cells, 
                                   Epetra_Vector& pressure_cells_dSdP, double dTp, Matrix_MFD* matrix);

  double ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot);
  void ComputePreconditionerMFD(const Epetra_Vector &u, Matrix_MFD* matrix,
                                double Tp, double dTp, bool flag_update_ML);

  void UpdateBoundaryConditions(double Tp, Epetra_Vector& p_faces);

  // linear problems and solvers
  void AssembleSteadyStateProblem_MFD(Matrix_MFD* matrix, bool add_preconditioner);
  void AssembleTransientProblem_MFD(Matrix_MFD* matrix, double dTp, Epetra_Vector& p, bool add_preconditioner);
  void SolveFullySaturatedProblem(double T, Epetra_Vector& u);
  void SolveTransientProblem(double T, double dT, Epetra_Vector& u);

  // io members
  void ProcessParameterList();
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringLinearSolver(const std::string name, int* max_itrs, double* tolerance);
  void ProcessStringRelativePermeability(const std::string name, int* method);
  void VerifyStringMualemBurdine(const std::string name);

  std::string FindStringPreconditioner(const Teuchos::ParameterList& list);
  std::string FindStringLinearSolver(const Teuchos::ParameterList& list);
  void ProcessSublistTimeIntegration(
    Teuchos::ParameterList& list, const std::string name,
    double* absolute_tol, double* relative_tol, double* residual_tol, int* max_itrs,
    double* T0, double* T1, double* dT0, double* dTmax);

  // water retention models
  void DerivedSdP(const Epetra_Vector& p, Epetra_Vector& dS);
  void DeriveSaturationFromPressure(const Epetra_Vector& p, Epetra_Vector& s);
  void DerivePressureFromSaturation(const Epetra_Vector& s, Epetra_Vector& p);

  // initization members
  void DeriveFaceValuesFromCellValues(const Epetra_Vector& ucells, Epetra_Vector& ufaces);
  void ClipHydrostaticPressure(const double pmin, Epetra_Vector& pressure_cells);
  void ClipHydrostaticPressure(const double pmin, const double s0, Epetra_Vector& pressure_cells);

  double CalculateRelaxationFactor(const Epetra_Vector& uold, const Epetra_Vector& unew);

  // control method
  void ResetErrorControl(double error) { error_control_ = error; }
  void ResetParameterList(const Teuchos::ParameterList& rp_list_new) { rp_list_ = rp_list_new; }
  void PrintStatistics() const;
  
  // access methods
  Teuchos::RCP<AmanziMesh::Mesh> mesh() { return mesh_; }
  const Epetra_Map& super_map() { return *super_map_; }
  AmanziGeometry::Point& gravity() { return gravity_; }

  // developement members
  void CalculateConsistentSaturation(const Epetra_Vector& flux, 
                                     const Epetra_Vector& ws_prev, Epetra_Vector& ws);

 public:
  int num_nonlinear_steps;

 private:
  Teuchos::ParameterList rp_list_;
  Teuchos::ParameterList solver_list_;
  Teuchos::ParameterList preconditioner_list_;

  AmanziGeometry::Point gravity_;
  double rho, mu;
  double atm_pressure;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map* super_map_;

  Teuchos::RCP<Epetra_Import> cell_importer_;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer_;

  AztecOO* solver;  // Linear solver data
  Matrix_MFD* matrix;
  Matrix_MFD* preconditioner;
  double convergence_tol;

  BDF2::Dae* bdf2_dae;  // Time intergrators
  BDF1Dae* bdf1_dae;
  int block_picard;
  int error_control_;

  int ti_method_sss;  // Parameters for steady-state solution
  std::string preconditioner_name_sss_;
  int num_itrs_sss, max_itrs_sss;
  double absolute_tol_sss, relative_tol_sss, residual_tol_sss;
  double T0_sss, T1_sss, dT0_sss, dTmax_sss;
  int initialize_with_darcy;

  int ti_method_trs;  // Parameters for transient solution
  std::string preconditioner_name_trs_;
  double absolute_tol_trs, relative_tol_trs, residual_tol_trs;
  int num_itrs_trs, max_itrs_trs;
  double T0_trs, T1_trs, dT0_trs, dTmax_trs;

  double absolute_tol, relative_tol;  // Generic parameters (sss or trs)
  int ti_method, num_itrs, max_itrs;

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.
  Teuchos::RCP<Epetra_Vector> rhs_faces;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM;

  BoundaryFunction* bc_pressure;  // Pressure BC.
  BoundaryFunction* bc_head;  // Static pressure head BC.
  BoundaryFunction* bc_flux;  // Outward mass flux BC.
  BoundaryFunction* bc_seepage;  // Seepage face BC.
  std::vector<int> bc_markers;  // Used faces are marked with boundary conditions.
  std::vector<double> bc_values;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  std::vector<AmanziGeometry::Point> Kgravity_unit;  // normalized vector Kg

  int Krel_method;  // method for calculating relative permeability
  Teuchos::RCP<Epetra_Vector> Krel_cells;  // realitive permeability 
  Teuchos::RCP<Epetra_Vector> Krel_faces;  // realitive permeability

  int mfd3d_method_, mfd3d_method_preconditioner_;
  bool is_matrix_symmetric;
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;

  double clip_saturation;  // initialization options

 private:
  void operator=(const Richards_PK& RPK);
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

