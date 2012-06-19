/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __DARCY_PK_HPP__
#define __DARCY_PK_HPP__

#include "Teuchos_RCP.hpp"

#include "Epetra_Vector.h"
#include "AztecOO.h"

#include "Mesh.hh"
#include "Point.hh"
#include "boundary_function.hh"
#include "domain_function.hh"
#include "tensor.hpp"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Darcy_PK : public Flow_PK {
 public:
  Darcy_PK(Teuchos::ParameterList& global_list, Teuchos::RCP<Flow_State> FS_MPC);
  ~Darcy_PK();

  // main methods
  void InitPK();
  void InitSteadyState(double T0, double dT0);
  void InitTransient(double T0, double dT0);
  void InitPicard(double T0) {};  // not used yet.

  double CalculateFlowDt() { return dT_desirable_; }
  int Advance(double dT); 
  int AdvanceToSteadyState();
  void InitializeAuxiliaryData();
  void InitializeSteadySaturated();

  void CommitState(Teuchos::RCP<Flow_State> FS);

  // methods required for time integration
  void fun(const double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs, double dT = 0.0) {};
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu) {};
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du) { return 0.0; }
  void update_precon(const double T, const Epetra_Vector& up, const double h, int& errc) {};
  void update_norm(double rtol, double atol) {};

  // other main methods
  void SetAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K);
  void AddTimeDerivativeSpecificStorage(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix_operator);
  void AddTimeDerivativeSpecificYield(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix_operator);
  void UpdateSpecificYield();

  // linear solvers
  void SolveFullySaturatedProblem(double T, Epetra_Vector& u);

  // io members
  void ProcessParameterList();
  void ProcessStringLinearSolver(const std::string name, int* max_itrs, double* tolerance);

  // control methods
  void PrintStatistics() const;
  void ResetParameterList(const Teuchos::ParameterList& dp_list_new) { dp_list_ = dp_list_new; }

  // access methods
  Epetra_Vector& ref_solution_faces() { return *solution_faces; }
  Epetra_Import& ref_face_importer() { return *face_importer_; }

  double rho() { return rho_; }
  double mu() { return mu_; }
  AmanziGeometry::Point& gravity() { return gravity_; }

 private:
  Teuchos::ParameterList dp_list_;
  Teuchos::ParameterList preconditioner_list_;
  Teuchos::ParameterList solver_list_;

  AmanziGeometry::Point gravity_;
  double rho_, mu_;
  double atm_pressure;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map* super_map_;
  int dim;

  Teuchos::RCP<Epetra_Import> cell_importer_;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer_;

  AztecOO* solver;
  Matrix_MFD* matrix_;
  Matrix_MFD* preconditioner_;

  int num_itrs_sss, max_itrs_sss;  // Parameters for steady state solution
  std::string preconditioner_name_sss_;
  double convergence_tol_sss, residual_sss;

  int num_itrs_trs;  // Parameters for transient solver
  double dT_desirable_;

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.
  Teuchos::RCP<Epetra_Vector> rhs_faces;

  BoundaryFunction* bc_pressure;  // Boundary conditions. 
  BoundaryFunction* bc_head;
  BoundaryFunction* bc_flux;
  BoundaryFunction* bc_seepage;
  std::vector<int> bc_markers;  // Used faces marked with boundary conditions
  std::vector<double> bc_values;

  DomainFunction* src_sink;  // Source and sink terms

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  Teuchos::RCP<Epetra_Vector> Krel_cells;  // realitive permeability 
  Teuchos::RCP<Epetra_Vector> Krel_faces;  // realitive permeability 

  int mfd3d_method;
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

