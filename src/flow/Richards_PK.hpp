/*
This is the flow component of the Amanzi code. 
License: BSD
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

#include "mfd3d.hpp"
#include "BDF2_Dae.hpp"

#include "Flow_State.hpp"
#include "Flow_PK.hpp"
#include "Flow_BC_Factory.hpp"
#include "Matrix_MFD.hpp"
#include "WaterRetentionModel.hpp"


namespace Amanzi {
namespace AmanziFlow {

class Interface_BDF2;  // forward declaration of class

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::RCP<Teuchos::ParameterList> rp_list_, Teuchos::RCP<Flow_State> FS_MPC);
  ~Richards_PK () { delete super_map_, solver, matrix, preconditioner, bc_pressure, bc_head, bc_flux; }

  // main methods
  void Init(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL);

  int advance(double dT); 
  int advance_to_steady_state();
  void commit_state() {};

  int advanceSteadyState_Picard();
  int advanceSteadyState_BackwardEuler();
  int advanceSteadyState_NOX();
 
  // other main methods
  void processParameterList();
  void setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K);
  void calculateRelativePermeability(const Epetra_Vector& p);
  void calculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p);
  void calculateRelativePermeabilityUpwindFlux(const Epetra_Vector& p, const Epetra_Vector& darcy_flux);

  void applyDiffusionMFD(Epetra_Vector& X, Epetra_Vector& Y);
  void applyAdvectionMFD(Epetra_Vector& X, Epetra_Vector& Y);

  void computePDot(const double T, const Epetra_Vector& p, Epetra_Vector& pdot);
  void computeFunctionalPTerm(const Epetra_Vector& p, Epetra_Vector& pdot, double T = 0.0);

  void computePreconditioner(const Epetra_Vector &X, const double h);

  void derivedSdP(const Epetra_Vector& p, Epetra_Vector& dS);
  void deriveVanGenuchtenSaturation(const Epetra_Vector& p, Epetra_Vector& s);
  void derivePressureFromSaturation(double s, Epetra_Vector& p);

  void addGravityFluxes_MFD(Matrix_MFD* matrix);
  void addTimeDerivative_MFD(Epetra_Vector& pressure_cells, Matrix_MFD* matrix);

  // control methods
  void print_statistics() const;
  
  // access methods
  Flow_State& get_FS() { return *FS; }
  Teuchos::RCP<AmanziMesh::Mesh> get_mesh() { return mesh_; }
  Matrix_MFD* get_matrix() { return matrix; }
  Epetra_Vector& get_solution_cells() { return *solution_cells; }
  Epetra_Vector& get_solution_faces() { return *solution_faces; }
  AmanziGeometry::Point& get_gravity() { return gravity; }
  double get_rho() { return rho; }
  double get_mu() { return mu; }

 private:
  Teuchos::RCP<Teuchos::ParameterList> rp_list;

  Teuchos::RCP<Flow_State> FS;
  AmanziGeometry::Point gravity;
  double rho, mu;
  double atm_pressure;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map *super_map_;
  int dim;

  Teuchos::RCP<Epetra_Import> cell_importer_;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer_;

  AztecOO* solver;
  Matrix_MFD* matrix;
  Matrix_MFD* preconditioner;

  int precon_freq;  // preconditioner update frequency
  int num_itrs, max_itrs;  // numbers of linear solver iterations
  double err_tol;  // errors in linear solver

  int method_sss;  // Method to reach the steady-state solution
  int num_itrs_sss, max_itrs_sss;  // number of non-linear iterations
  double err_tol_sss, residual_sss;

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.
  Teuchos::RCP<Epetra_Vector> rhs_faces;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM;
  Interface_BDF2* solver_BDF;
  BDF2::Dae* time_stepper;

  BoundaryFunction *bc_pressure;  // Pressure Dirichlet b.c., excluding static head
  BoundaryFunction *bc_head;  // Static pressure head b.c.; also Dirichlet-type
  BoundaryFunction *bc_flux;  // Outward mass flux b.c.
  std::vector<int> bc_markers;  // Used faces marked with boundary conditions
  std::vector<double> bc_values;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  Teuchos::RCP<Epetra_Vector> Krel_cells;  // realitive permeability 
  Teuchos::RCP<Epetra_Vector> Krel_faces;  // realitive permeability 

  bool upwind_Krel;
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

