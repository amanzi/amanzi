#ifndef __RICHARDSPROBLEM_HH__
#define __RICHARDSPROBLEM_HH__

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

namespace Amanzi {

class BoundaryFunction; // forward declaration

class RichardsProblem
{
public:
  RichardsProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                  Teuchos::ParameterList& richards_plist);
  ~RichardsProblem();

  // initialization
  void InitializeProblem(Teuchos::ParameterList& plist);

  // set independent variables
  void SetPermeability(double k);
  void SetPermeability(const Epetra_Vector &k);
  void SetPorosity(double k);
  void SetPorosity(const Epetra_Vector &k);
  void SetFluidDensity(double rho);
  void SetFluidViscosity(double mu);
  void SetGravity(const double g[3]);
  void SetGravity(double g); // note this is g in the negative-z-direction

  // derive things
  void DeriveSaturation(const Epetra_Vector& P, Epetra_Vector& S) const;
  void DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l1_error) const;
  void DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const;

  // timestepping
  void Compute_udot(const double t, const Epetra_Vector& u, Epetra_Vector &udot);
  void SetInitialPressureProfileCells(double height, Teuchos::RCP<Epetra_Vector> &pressure);
  void SetInitialPressureProfileFaces(double height, Teuchos::RCP<Epetra_Vector> &pressure);
  void SetInitialPressureProfileFromSaturationCells(double saturation,
                                                    Teuchos::RCP<Epetra_Vector> &pressure);
  void SetInitialPressureProfileFromSaturationFaces(double saturation,
                                                    Teuchos::RCP<Epetra_Vector> &pressure);
  // access
  double GetFluidDensity() const { return rho_; }
  double GetFluidViscosity() const { return mu_; }
  Teuchos::RCP<double*> GetGravity() { return gvec_; }

  Teuchos::RCP<Epetra_Vector> GetPorosity() { return phi_; }
  Teuchos::RCP<Epetra_Vector> GetCellVolumes() { return cell_volumes_; }

  Epetra_Vector* CreateCellView(const Epetra_Vector&) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector&) const;

  const Epetra_Map& Map() const { return *dof_map_; }
  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }
  DiffusionMatrix& Matrix() const { return *D_; }

  // residuals/preconditioners for model evaluator
  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time);
  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F);

  void ComputePrecon(const Epetra_Vector &X);
  void ComputePrecon(const Epetra_Vector &X, double h);
  void dSofP(const Epetra_Vector& P, Epetra_Vector& dS);
  Epetra_Operator& Precon() { return *precon_; }
private:

  // START BELOW SHOULD BE PRIVATE??
  // intermediate steps
  void ComputeRelPerm(const Epetra_Vector&, Epetra_Vector&) const;
  void ComputeUpwindRelPerm(const Epetra_Vector& Pcell, const Epetra_Vector& Pface,
                            Epetra_Vector& k_rel) const;
  void UpdateVanGenuchtenRelativePermeability(const Epetra_Vector& P);

  // END BELOW SHOULD BE PRIVATE??

  // access

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }

  Teuchos::RCP<Epetra_Map> create_dof_map_(const Epetra_Map&, const Epetra_Map&) const;
  void validate_boundary_conditions() const;
  DiffusionMatrix* create_diff_matrix_(Teuchos::RCP<AmanziMesh::Mesh>&) const;
  void init_mimetic_disc_(Teuchos::RCP<AmanziMesh::Mesh>&, std::vector<MimeticHexLocal>&) const;
  void upwind_rel_perm_(const Epetra_Vector&, Epetra_Vector&);

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Epetra_Map> dof_map_;
  Teuchos::RCP<Epetra_Import> face_importer_;
  Teuchos::RCP<Epetra_Import> cell_importer_;

  std::vector<double> k_; // spatially variable permeability
  std::vector<double> k_rl_;  // relative permeability
  double p_atm_;    // atmospheric pressure
  double rho_;      // fluid density
  double mu_;       // fluid viscosity
  double gravity_;  // gravitational acceleration (positive coef, directed in -z direction)
  Teuchos::RCP<double*> gvec_;  // set to {0,0,-gravity_}

  Teuchos::RCP<Epetra_Vector> phi_;
  Teuchos::RCP<Epetra_Vector> cell_volumes_;

  Teuchos::RCP<BoundaryFunction> bc_press_;  // Pressure Dirichlet conditions, excluding static head
  Teuchos::RCP<BoundaryFunction> bc_head_;   // Static pressure head conditions; also Dirichlet-type
  Teuchos::RCP<BoundaryFunction> bc_flux_;   // Outward mass flux conditions

  std::vector<MimeticHexLocal>  MD_;
  MimeticHex *md_;
  bool upwind_k_rel_;

  DiffusionPrecon *precon_;

  Teuchos::RCP<DiffusionMatrix> D_;


  std::vector< Teuchos::RCP<WaterRetentionBaseModel> > WRM_;

};

} // close namespace Amanzi

#endif
