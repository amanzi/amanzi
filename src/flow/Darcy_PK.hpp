/*
This is the flow component of the Amanzi code. 
License: BSD
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
#include "boundary-function.hh"
#include "tensor.hpp"

#include "Flow_PK.hpp"
#include "Flow_State.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Darcy_PK : public Flow_PK {
 public:
  Darcy_PK(Teuchos::ParameterList& flow_list, Teuchos::RCP<Flow_State> FS_MPC);
  ~Darcy_PK();

  // main methods
  void InitPK(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL);
  void InitSteadyState() {};  // (lipnikov@lanl.gov)
  void InitTransient() {};

  int advance(double dT) { return 0; } 
  int advance_to_steady_state();
  void commit_state(Teuchos::RCP<Flow_State> FS) {};

  // required methods
  void fun(const double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs) {};
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu) {};
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du) { return 0.0; }
  void update_precon(const double T, const Epetra_Vector& up, const double h, int& errc) {};

  // other main methods
  void processParameterList();
  void setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K);
  void deriveDarcyVelocity(Epetra_MultiVector& velocity);

  // control methods
  void print_statistics() const;
  void resetParameterList(const Teuchos::ParameterList& dp_list_new) { dp_list = dp_list_new; }

  // access methods
  Epetra_Vector& ref_solution_faces() { return *solution_faces; }

  double rho() { return rho_; }
  double mu() { return mu_; }
  AmanziGeometry::Point& gravity() { return gravity_; }

 private:
  Teuchos::ParameterList dp_list;

  AmanziGeometry::Point gravity_;
  double rho_, mu_;
  double atm_pressure;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map* super_map_;
  int dim;

  Teuchos::RCP<Epetra_Import> cell_importer_;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer_;

  AztecOO* solver;
  Matrix_MFD* matrix;
  Matrix_MFD* preconditioner;

  int num_itrs_sss, max_itrs_sss;  // Parameters for steady state solution
  double convergence_tol_sss, residual_sss;

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.
  Teuchos::RCP<Epetra_Vector> rhs_faces;

  BoundaryFunction* bc_pressure;  // Pressure Dirichlet b.c., excluding static head
  BoundaryFunction* bc_head;  // Static pressure head b.c.; also Dirichlet-type
  BoundaryFunction* bc_flux;  // Outward mass flux b.c.
  std::vector<int> bc_markers;  // Used faces marked with boundary conditions
  std::vector<double> bc_values;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  Teuchos::RCP<Epetra_Vector> Krel_cells;  // realitive permeability 
  Teuchos::RCP<Epetra_Vector> Krel_faces;  // realitive permeability 

  int mfd3d_method;
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

