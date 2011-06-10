#ifndef __RICHARDSPROBLEM_H__
#define __RICHARDSPROBLEM_H__

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh_maps_base.hh"

#include "FlowBC.hpp"
#include "DiffusionMatrix.hpp"
#include "DiffusionPrecon.hpp"
#include "MimeticHexLocal.hpp"
#include "MimeticHex.hpp"
#include "WaterRetentionBaseModel.hpp"
#include "Flow_State.hpp"

class RichardsProblem 
{
public:

  RichardsProblem(const Teuchos::RCP<Mesh_maps_base> &mesh,
	       Teuchos::ParameterList&, const Teuchos::RCP<FlowBC> &bc);
  ~RichardsProblem();


  // Set the constant value of fluid density.
  void SetFluidDensity(double rho);

  // Set the constant value of fluid viscosity.
  void SetFluidViscosity(double mu);

  // Set the gravitational acceleration (3-vector).
  void SetGravity(const double g[3]);

  // Sets a constant (scalar) permeability.
  void SetPermeability(double k);

  // Sets a spatially variable (scalar) permeability, one value per cell.
  //void SetPermeability(const std::vector<double> &k);
  void SetPermeability(const Epetra_Vector &k);

  void SetFlowState( Teuchos::RCP<const Flow_State> FS_ );

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
  void GetGravity(double g[]) const { for(int i = 0; i < 3; ++i) g[i] = g_[i]; }

  void Compute_udot(const double t, const Epetra_Vector& u, Epetra_Vector &udot);
  
  const Epetra_Vector* cell_vols() { return cell_volumes; }

  void SetInitialPressureProfileCells(double height, Epetra_Vector *pressure);
  void SetInitialPressureProfileFaces(double height, Epetra_Vector *pressure);
  void SetInitialPressureProfileFromSaturationCells(double saturation, Epetra_Vector *pressure);
  void SetInitialPressureProfileFromSaturationFaces(double saturation, Epetra_Vector *pressure);

private:

  Teuchos::RCP<Mesh_maps_base> mesh_;
  Teuchos::RCP<FlowBC> bc_;
  Epetra_Map *dof_map_;
  Epetra_Import *face_importer_;
  Epetra_Import *cell_importer_;

  double rho_;  // constant fluid density
  double mu_;   // constant fluid viscosity
  double g_[3]; // gravitational acceleration
  std::vector<double> k_; // spatially variable permeability
  std::vector<double> k_rl_;  // relative permeability
  //double vG_m_;     // van Genuchten m
  //double vG_n_;     // van Genuchten n = 1/(1-vG_m_)
  //double vG_alpha_; // van Genuchten alpha
  //double vG_sr_;    // van Genuchten effective saturation
  double p_atm_;    // atmospheric pressure
  
  std::vector<MimeticHexLocal>  MD;
  MimeticHex *md_;

  DiffusionPrecon *precon_;

  Teuchos::RCP<DiffusionMatrix> D_;

  Epetra_Vector* cell_volumes;
  
  std::vector<Teuchos::RCP<WaterRetentionBaseModel> > WRM;

  Teuchos::RCP<const Flow_State> FS;  

private:  // Auxillary functions

  Epetra_Map* create_dof_map_(const Epetra_Map&, const Epetra_Map&) const;
  DiffusionMatrix* create_diff_matrix_(const Teuchos::RCP<Mesh_maps_base>&, const Teuchos::RCP<FlowBC>&) const;
  void init_mimetic_disc_(Mesh_maps_base&, std::vector<MimeticHexLocal>&) const;
  void apply_BC_initial_(Epetra_Vector&, double);
  void apply_BC_final_(Epetra_Vector&);
  void face_centroid_(int, double[]);


};

#endif
