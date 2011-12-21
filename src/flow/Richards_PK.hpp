/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __RICHARDS_PK_HPP__
#define __RICHARDS_PK_HPP__

#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "mfd3d.hpp"
#include "BDF2_Dae.hpp"

#include "Flow_State.hpp"
#include "Flow_PK.hpp"
#include "Flow_BC_Factory.hpp"
#include "Matrix_MFD.hpp"
#include "WaterRetentionModel.hpp"
#include "Interface_BDF2.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::ParameterList& rp_list_, Teuchos::RCP<Flow_State> FS_MPC);
  ~Richards_PK () { delete super_map_, solver, matrix, preconditioner, bc_pressure, bc_head, bc_flux; }

  // main methods
  int advance(double dT); 
  int advance_to_steady_state();
  void commit_state(Teuchos::RCP<Flow_State>) {};
  int init(double t0, double h0);

  // other main methods
  void process_parameter_list();
  void populate_absolute_permeability_tensor(std::vector<WhetStone::Tensor>& K);
  void calculate_relative_permeability(std::vector<double>& Krel);

  void compute_pdot(const double T, const Epetra_Vector& p, Epetra_Vector &pdot);

  // control methods
  void print_statistics() const;
  
 private:
  void Init(Matrix_MFD* matrix_ = NULL, Matrix_MFD* preconditioner_ = NULL);
  Teuchos::ParameterList rp_list;

  Teuchos::RCP<Flow_State> FS;
  AmanziGeometry::Point gravity;
  double rho, mu;
  double atm_pressure;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map *super_map_;
  int dim;

  Teuchos::RCP<Epetra_Import> cell_importer_;  // parallel communicators
  Teuchos::RCP<Epetra_Import> face_importer_;

  Matrix_MFD* matrix;
  Matrix_MFD* preconditioner;

  int precon_freq;  // preconditioner update frequency
  int num_itrs, max_itrs;  // numbers of linear solver iterations
  double err_tol, residual;  // errors in linear solver

  Teuchos::RCP<Epetra_Vector> solution;  // global solution
  Teuchos::RCP<Epetra_Vector> solution_cells;  // cell-based pressures
  Teuchos::RCP<Epetra_Vector> solution_faces;  // face-base pressures
  Teuchos::RCP<Epetra_Vector> rhs;  // It has same size as solution.

  std::vector<Teuchos::RCP<WaterRetentionModel> > WRM;
  Interface_BDF2* solver;
  BDF2::Dae* time_stepper;

  BoundaryFunction *bc_pressure;  // Pressure Dirichlet b.c., excluding static head
  BoundaryFunction *bc_head;  // Static pressure head b.c.; also Dirichlet-type
  BoundaryFunction *bc_flux;  // Outward mass flux b.c.
  std::vector<int> bc_markers;  // Used faces marked with boundary conditions
  std::vector<double> bc_values;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability
  std::vector<double> Krel;  // realitive permeability
  bool upwind_Krel;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

/*
  void ComputeRelPerm(const Epetra_Vector&, Epetra_Vector&) const;
  void ComputeUpwindRelPerm(const Epetra_Vector&, const Epetra_Vector&, Epetra_Vector&) const;
  void UpdateVanGenuchtenRelativePermeability(const Epetra_Vector &P);
  void DeriveVanGenuchtenSaturation(const Epetra_Vector &P, Epetra_Vector &S);
  void dSofP(const Epetra_Vector &P, Epetra_Vector &dS);
  
  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time = 0.0);
  void ComputePrecon(const Epetra_Vector &X);
  void ComputePrecon(const Epetra_Vector &X, const double h);

  void compute_udot(const double t, const Epetra_Vector& u, Epetra_Vector& udot);
  
  DiffusionMatrix& Matrix() const { return *D_; }
  void DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l1_error) const;
  void DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const;

  std::vector<Teuchos::RCP<WaterRetentionBaseModel> > WRM;
  void upwind_rel_perm_(const Epetra_Vector&, Epetra_Vector&);
*/
