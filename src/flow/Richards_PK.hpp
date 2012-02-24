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
  Richards_PK(Teuchos::ParameterList& flow_list, Teuchos::RCP<Flow_State> FS_MPC);
  ~Richards_PK();

  // main methods
  void InitPK(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL);
  void InitSteadyState(double T0, double dT0);
  void InitTransient(double T0, double dT0);

  int advance(double dT_MPC); 
  int advanceToSteadyState();
  void commitState(Teuchos::RCP<Flow_State> FS);
  void commitStateForTransport(Teuchos::RCP<Flow_State> FS);
  void deriveDarcyVelocity(const Epetra_Vector& flux, Epetra_MultiVector& velocity);

  int advanceSteadyState_Picard();
  int advanceSteadyState_BackwardEuler();
  int advanceSteadyState_ForwardEuler();
  int advanceSteadyState_BDF2();
 
  // methods required for time integration
  void fun(double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& rhs, double dT = 0.0);
  void precon(const Epetra_Vector& u, Epetra_Vector& Hu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_precon(double T, const Epetra_Vector& u, double dT, int& ierr);

  // other main methods
  void processParameterList();
  void setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K);
  void calculateRelativePermeabilityCell(const Epetra_Vector& p);
  void calculateRelativePermeabilityFace(const Epetra_Vector& p);
  void calculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p);
  void calculateRelativePermeabilityUpwindFlux(const Epetra_Vector& p, const Epetra_Vector& flux);
  void calculateRelativePermeabilityArithmeticMean(const Epetra_Vector& p);

  void addTimeDerivative_MFD(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix);
  void addTimeDerivative_MFDfake(Epetra_Vector& pressure_cells, double dTp, Matrix_MFD* matrix);

  double computeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot);
  void computePreconditionerMFD(
      const Epetra_Vector &u, Matrix_MFD* matrix, double Tp, double dTp, bool flag_update_ML);
  double errorSolutionDiff(const Epetra_Vector& uold, const Epetra_Vector& unew);

  void derivedSdP(const Epetra_Vector& p, Epetra_Vector& dS);
  void deriveSaturationFromPressure(const Epetra_Vector& p, Epetra_Vector& s);
  void derivePressureFromSaturation(const Epetra_Vector& s, Epetra_Vector& p);

  void deriveFaceValuesFromCellValues(const Epetra_Vector& ucells, Epetra_Vector& ufaces);

  // control methods
  void resetParameterList(const Teuchos::ParameterList& rp_list_new) { rp_list = rp_list_new; }
  void print_statistics() const;
  
  // access methods
  Teuchos::RCP<AmanziMesh::Mesh> mesh() { return mesh_; }
  const Epetra_Map& super_map() { return *super_map_; }
  AmanziGeometry::Point& gravity() { return gravity_; }

 public:
  int num_nonlinear_steps;

 private:
  Teuchos::ParameterList rp_list;

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
  int max_itrs;
  double convergence_tol;

  int method_sss;  // Parameters for steady-state solution
  int num_itrs_sss, max_itrs_sss;
  double absolute_tol_sss, relative_tol_sss, convergence_tol_sss;
  double T0_sss, T1_sss, dT0_sss, dTmax_sss;

  BDF2::Dae* bdf2_dae;  // Parameters for transient solution
  int method_trs;
  double absolute_tol_trs, relative_tol_trs, convergence_tol_trs;
  int num_itrs_trs, max_itrs_trs;
  double T0_trs, T1_trs, dT0_trs, dTmax_trs;

  double absolute_tol, relative_tol;  // Generic parameters (sss or trs)

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.
  Teuchos::RCP<Epetra_Vector> rhs_faces;

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM;

  BoundaryFunction* bc_pressure;  // Pressure Dirichlet b.c., excluding static head
  BoundaryFunction* bc_head;  // Static pressure head b.c.; also Dirichlet-type
  BoundaryFunction* bc_flux;  // Outward mass flux b.c.
  std::vector<int> bc_markers;  // Used faces are marked with boundary conditions.
  std::vector<double> bc_values;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  Teuchos::RCP<Epetra_Vector> Krel_cells;  // realitive permeability 
  Teuchos::RCP<Epetra_Vector> Krel_faces;  // realitive permeability 

  int Krel_method;  // method for calculating relative permeability
  int mfd3d_method;
  bool is_matrix_symmetric;
  Teuchos::RCP<Epetra_IntVector> upwind_cell, downwind_cell;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

