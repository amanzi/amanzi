/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __RICHARDS_PK_HPP__
#define __RICHARDS_PK_HPP__

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Flow_State.hpp"
#include "Flow_PK.hpp"
#include "Interface_BDF2.hpp"

namespace Amanzi {
namespace AmanziFlow {

class Richards_PK : public Flow_PK {
 public:
  Richards_PK(Teuchos::ParameterList&, const Teuchos::RCP<const Flow_State>);
  ~Richards_PK ();

  // major memebr
  int init(double t0, double h0);
  int advance(double dT);
  int advance_to_steady_state();
  void commit_state(Teuchos::RCP<Flow_State>) {};  // pointer to state is know

  // access members
  const Epetra_Vector& get_pressure() const { return *pressure_cells; }
  const Epetra_Vector& get_flux() const { return *richards_flux; }
  void get_velocity(Epetra_MultiVector &q) const { problem->DeriveDarcyVelocity(*solution, q); }
  void get_saturation(Epetra_Vector &s) const;
  double get_flow_dT() { return hnext; }
  
 private:
  Teuchos::RCP<const Flow_State> FS;
  Teuchos::ParameterList &richards_plist;
 
  RichardsModelEvaluator *RME;
  
  BDF2::Dae *time_stepper;

  Epetra_Vector *solution;   // full cell/face solution
  Epetra_Vector *pressure_cells;   // cell pressures
  Epetra_Vector *pressure_faces;
  Epetra_Vector *richards_flux; // Darcy face fluxes

  int max_itr;      // max number of linear solver iterations
  double err_tol;   // linear solver convergence error tolerance
  int precon_freq;  // preconditioner update frequency

  double ss_t0, ss_t1, ss_h0, ss_z;

  double h, hnext;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

/*
class RichardsProblem {
 public:
  RichardsProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh, 
                  const Teuchos::ParameterList& parameter_list);
  ~RichardsProblem();

  // explicit initialization
  void set_absolute_permeability(double k);
  void set_absolute_permeability(const Epetra_Vector &k);

  void set_pressure_cells(double height, Epetra_Vector *pressure);
  void set_pressure_faces(double height, Epetra_Vector *pressure);
  void set_initial_pressure_from_saturation_cells(double saturation, Epetra_Vector *pressure);
  void set_initial_pressure_from_saturation_faces(double saturation, Epetra_Vector *pressure);

  void set_flow_state(Teuchos::RCP<const Flow_State> FS_) { FS = FS_; } 

  void ComputeRelPerm(const Epetra_Vector&, Epetra_Vector&) const;
  void ComputeUpwindRelPerm(const Epetra_Vector&, const Epetra_Vector&, Epetra_Vector&) const;
  void UpdateVanGenuchtenRelativePermeability(const Epetra_Vector &P);
  void DeriveVanGenuchtenSaturation(const Epetra_Vector &P, Epetra_Vector &S);
  void dSofP(const Epetra_Vector &P, Epetra_Vector &dS);
  
  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time = 0.0);

  void ComputePrecon(const Epetra_Vector &X);
  void ComputePrecon(const Epetra_Vector &X, const double h);

  void compute_udot(const double t, const Epetra_Vector& u, Epetra_Vector& udot);
  
  // maps
  const Epetra_Map& Map() const { return *dof_map_; }
  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }

  Epetra_Vector* CreateCellView(const Epetra_Vector&) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector&) const;

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }

  DiffusionMatrix& Matrix() const { return *D_; }
  void DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l1_error) const;
  void DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const;

  // access
  void GetFluidDensity(double &rho) const { rho = rho_; }
  void GetFluidViscosity(double &mu) const { mu = mu_; }
  void GetGravity(double g[]) const { for(int i=0; i<3; ++i) g[i] = gvec_[i]; }

  Epetra_Operator& Precon() const { return *precon_; }

  const Epetra_Vector* cell_vols() { return cell_volumes; }

 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map *dof_map_;
  Epetra_Import *face_importer_;
  Epetra_Import *cell_importer_;

  std::vector<double> k_;  // spatially variable permeability
  std::vector<double> k_rl_;  // relative permeability
  double p_atm_;  // atmospheric pressure
  double rho_;  // constant fluid density
  double mu_;  // constant fluid viscosity
  double gravity_;  // gravitational acceleration (positive coef, directed in -z direction)
  double gvec_[3];  // set to {0,0,-gravity_}
  
  Amanzi::BoundaryFunction *bc_press_;  // Pressure Dirichlet conditions, excluding static head
  Amanzi::BoundaryFunction *bc_head_;   // Static pressure head conditions; also Dirichlet-type
  Amanzi::BoundaryFunction *bc_flux_;   // Outward mass flux conditions
  
  std::vector<MimeticHexLocal>  MD;
  MimeticHex *md_;
  bool upwind_k_rel_;

  DiffusionPrecon *precon_;

  Teuchos::RCP<DiffusionMatrix> D_;

  std::vector<Teuchos::RCP<WaterRetentionBaseModel> > WRM;

  Teuchos::RCP<const Flow_State> FS;  

 private:
  Epetra_Map* create_dof_map_(const Epetra_Map&, const Epetra_Map&) const;
  void validate_boundary_conditions() const;
  DiffusionMatrix* create_diff_matrix_(const Teuchos::RCP<AmanziMesh::Mesh>&) const;
  void init_mimetic_disc_(AmanziMesh::Mesh&, std::vector<MimeticHexLocal>&) const;
  void upwind_rel_perm_(const Epetra_Vector&, Epetra_Vector&);
};

*/
