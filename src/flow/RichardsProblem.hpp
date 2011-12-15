#ifndef __RICHARDSPROBLEM_H__
#define __RICHARDSPROBLEM_H__

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"

#include "DiffusionMatrix.hpp"
#include "DiffusionPrecon.hpp"
#include "MimeticHexLocal.hpp"
#include "MimeticHex.hpp"
#include "WaterRetentionBaseModel.hpp"
#include "Flow_State.hpp"

namespace Amanzi {

class BoundaryFunction; // forward declaration

class RichardsProblem {
 public:
  RichardsProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh, 
                  const Teuchos::ParameterList& parameter_list);
  ~RichardsProblem();

  // Sets a constant (scalar) permeability.
  void SetPermeability(double k);

  // Sets a spatially variable (scalar) permeability, one value per cell.
  void SetPermeability(const Epetra_Vector &k);

  void SetFlowState(Teuchos::RCP<const Flow_State> FS_);

  void ComputeRelPerm(const Epetra_Vector&, Epetra_Vector&) const;
  void ComputeUpwindRelPerm(const Epetra_Vector&, const Epetra_Vector&, Epetra_Vector&) const;
  void UpdateVanGenuchtenRelativePermeability(const Epetra_Vector &P);
  void DeriveVanGenuchtenSaturation(const Epetra_Vector &P, Epetra_Vector &S);
  void dSofP(const Epetra_Vector &P, Epetra_Vector &dS);
  
  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time);
  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F);

  void ComputePrecon(const Epetra_Vector &X);
  void ComputePrecon(const Epetra_Vector &X, const double h);

  Epetra_Operator& Precon() const { return *precon_; }

  const Epetra_Map& Map() const { return *dof_map_; }

  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }

  Epetra_Vector* CreateCellView(const Epetra_Vector&) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector&) const;

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }

  DiffusionMatrix& Matrix() const { return *D_; }

  void DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l1_error) const;

  void DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const;

  void GetFluidDensity(double &rho) const { rho = rho_; }
  void GetFluidViscosity(double &mu) const { mu = mu_; }
  void GetGravity(double g[]) const { for(int i = 0; i < 3; ++i) g[i] = gvec_[i]; }

  void Compute_udot(const double t, const Epetra_Vector& u, Epetra_Vector &udot);
  
  const Epetra_Vector* cell_vols() { return cell_volumes; }

  void SetInitialPressureProfileCells(double height, Epetra_Vector *pressure);
  void SetInitialPressureProfileFaces(double height, Epetra_Vector *pressure);
  void SetInitialPressureProfileFromSaturationCells(double saturation, Epetra_Vector *pressure);
  void SetInitialPressureProfileFromSaturationFaces(double saturation, Epetra_Vector *pressure);

 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Epetra_Map *dof_map_;
  Epetra_Import *face_importer_;
  Epetra_Import *cell_importer_;

  std::vector<double> k_; // spatially variable permeability
  std::vector<double> k_rl_;  // relative permeability
  double p_atm_;    // atmospheric pressure
  double rho_;      // fluid density
  double mu_;       // fluid viscosity
  double gravity_;  // gravitational acceleration (positive coef, directed in -z direction)
  double gvec_[3];  // set to {0,0,-gravity_}
  
  BoundaryFunction *bc_press_;  // Pressure Dirichlet conditions, excluding static head
  BoundaryFunction *bc_head_;   // Static pressure head conditions; also Dirichlet-type
  BoundaryFunction *bc_flux_;   // Outward mass flux conditions
  
  std::vector<MimeticHexLocal>  MD;
  MimeticHex *md_;
  bool upwind_k_rel_;

  DiffusionPrecon *precon_;

  Teuchos::RCP<DiffusionMatrix> D_;

  Epetra_Vector* cell_volumes;
  
  std::vector<Teuchos::RCP<WaterRetentionBaseModel> > WRM;

  Teuchos::RCP<const Flow_State> FS;  

 private:  // Auxillary functions
  Epetra_Map* create_dof_map_(const Epetra_Map&, const Epetra_Map&) const;
  void validate_boundary_conditions() const;
  DiffusionMatrix* create_diff_matrix_(const Teuchos::RCP<AmanziMesh::Mesh>&) const;
  void init_mimetic_disc_(AmanziMesh::Mesh&, std::vector<MimeticHexLocal>&) const;
  void upwind_rel_perm_(const Epetra_Vector&, Epetra_Vector&);

};

} // close namespace Amanzi

#endif
