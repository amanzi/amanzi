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

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::ParameterList& rp_list_, Teuchos::RCP<Flow_State> FS_MPC);
  ~Richards_PK ();

  // main methods
  void Init(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL);

  int advance(double dT); 
  int advance_to_steady_state();
  void commit_state(Teuchos::RCP<Flow_State> FS) {};

  int advanceSteadyState_Picard();
  int advanceSteadyState_BackwardEuler();
  int advanceSteadyState_ForwardEuler();
  int advanceSteadyState_BDF2();
 
  // required BDF2 methods
  void fun(const double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs);
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_precon(const double T, const Epetra_Vector& u, const double dT, int& ierr);

  // other main methods
  void processParameterList();
  void setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K);
  void calculateRelativePermeability(const Epetra_Vector& p);
  void calculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p);
  void calculateRelativePermeabilityUpwindFlux(const Epetra_Vector& p, const Epetra_Vector& darcy_flux);

  void applyDiffusionMFD(Epetra_Vector& X, Epetra_Vector& Y);
  void applyAdvectionMFD(Epetra_Vector& X, Epetra_Vector& Y);

  void addTimeDerivative_MFD(Epetra_Vector& pressure_cells, Matrix_MFD* matrix);

  double computeUDot(const double T, const Epetra_Vector& u, Epetra_Vector& udot);
  void computePreconditionerMFD(const Epetra_Vector &u, Matrix_MFD* matrix);

  void derivedSdP(const Epetra_Vector& p, Epetra_Vector& dS);
  void deriveVanGenuchtenSaturation(const Epetra_Vector& p, Epetra_Vector& s);
  void derivePressureFromSaturation(double s, Epetra_Vector& p);

  // control methods
  inline bool set_standalone_mode(bool mode) { standalone_mode = mode; } 
  void resetParameterList(const Teuchos::ParameterList& rp_list_new) { rp_list = rp_list_new; }
  void print_statistics() const;
  
  // access methods
  Teuchos::RCP<AmanziMesh::Mesh> get_mesh() { return mesh_; }
  const Epetra_Map& get_super_map() { return *super_map_; }

  Flow_State& get_FS() { return *FS; }
  Matrix_MFD* get_matrix() { return matrix; }

  Epetra_Vector& get_solution() { return *solution; }
  Epetra_Vector& get_solution_cells() { return *solution_cells; }
  Epetra_Vector& get_solution_faces() { return *solution_faces; }
  AmanziGeometry::Point& get_gravity() { return gravity; }

  double get_rho() { return rho; }
  double get_mu() { return mu; }

  std::vector<int>& get_bc_markers() { return bc_markers; }
  std::vector<double>& get_bc_values() { return bc_values; }
  Epetra_Vector& get_Krel_cells() { return *Krel_cells; }

 private:
  Teuchos::ParameterList rp_list;

  Teuchos::RCP<Flow_State> FS;
  AmanziGeometry::Point gravity;
  double rho, mu;
  double atm_pressure;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map* super_map_;
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

  BDF2::Dae* bdf2_dae;  // BDF2 time integration method
  double absolute_tol_bdf, relative_tol_bdf;
  int itrs_bdf2;

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.
  Teuchos::RCP<Epetra_Vector> rhs_faces;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM;

  BoundaryFunction *bc_pressure;  // Pressure Dirichlet b.c., excluding static head
  BoundaryFunction *bc_head;  // Static pressure head b.c.; also Dirichlet-type
  BoundaryFunction *bc_flux;  // Outward mass flux b.c.
  std::vector<int> bc_markers;  // Used faces marked with boundary conditions
  std::vector<double> bc_values;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  Teuchos::RCP<Epetra_Vector> Krel_cells;  // realitive permeability 
  Teuchos::RCP<Epetra_Vector> Krel_faces;  // realitive permeability 

  bool flag_upwind;
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;

  int verbosity;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

