/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __PERMAFROSTPROBLEM_HH__
#define __PERMAFROSTPROBLEM_HH__

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"

#include "DiffusionMatrix.hpp"
#include "DiffusionPrecon.hpp"
#include "MimeticHexLocal.hpp"
#include "MimeticHex.hpp"
#include "SaturationCurve.hh"
#include "SaturationCurveFactory.hh"

namespace Amanzi {

class BoundaryFunction; // forward declaration

class RichardsProblem {
public:
  RichardsProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                  Teuchos::ParameterList& richards_plist);
  ~RichardsProblem();

  // initialization
  void InitializeProblem(Teuchos::ParameterList& plist);

  // set independent variables
  // NOTE: some of these likely should never get used, as they should get calculated!
  void SetPorosity(const Epetra_Vector &phi) { *phi_ = phi; }
  void SetTemperature(const Epetra_Vector &temp) { *temp_ = temp; }

  void SetGasPressure(const Epetra_Vector &pressure_gas) { *pressure_gas_ = pressure_gas; }
  void SetGasSaturation(const Epetra_Vector &sat_gas) { *sat_gas_ = sat_gas; }
  void SetLiquidSaturation(const Epetra_Vector &sat_liquid) { *sat_liquid_ = sat_liquid; }
  void SetIceSaturation(const Epetra_Vector &sat_ice) { *sat_ice_ = sat_ice; }

  void SetPermeability(const Epetra_Vector &k) { *k_ = k; }
  void SetLiquidRelativePermeability(const Epetra_Vector &k) { *k_rel_liquid_ = k; }

  void SetGasDensity(const Epetra_Vector &rho_gas) { *rho_gas_ = rho_gas; }
  void SetLiquidDensity(const Epetra_Vector &rho_liquid) { *rho_liquid_ = rho_liquid; }
  void SetIceDensity(const Epetra_Vector &rho_ice) { *rho_ice_ = rho_ice; }
  void SetLiquidViscosity(const Epetra_Vector &mu_liquid) { *mu_liquid_ = mu_liquid; }

  void SetGravity(const double g[3]);
  void SetGravity(double g); // note this is g in the negative-z-direction

  // derive things
  // - mesh based
  void DeriveCellHeights(Teuchos::RCP<Epetra_Vector> &cell_heights) const;

  // - intermediate steps
  void DerivePermeability(const Epetra_Vector &phi,
                          Teuchos::RCP<Epetra_Vector> &perm) const;
  void DeriveLiquidRelPerm(const Epetra_Vector &saturation,
                           Teuchos::RCP<Epetra_Vector> &rel_perm) const;
  void DeriveGasDensity(const Epetra_Vector &pressure, const Epetra_Vector &temp,
                        Teuchos::RCP<Epetra_Vector> &rho) const;
  void DeriveLiquidDensity(const Epetra_Vector &pressure, const Epetra_Vector &temp,
                           Teuchos::RCP<Epetra_Vector> &rho) const;
  void DeriveIceDensity(const Epetra_Vector &pressure, const Epetra_Vector &temp,
                           Teuchos::RCP<Epetra_Vector> &rho) const;
  void DeriveLiquidViscosity(const Epetra_Vector &pressure, const Epetra_Vector &temp,
                             Teuchos::RCP<Epetra_Vector> &viscosity) const;
  void DeriveSaturation(const Epetra_Vector &pressure_gas,
                        const Epetra_Vector &pressure_liquid,
                        const Epetra_Vector &temp,
                        const Epetra_Vector &density_ice,
                        Teuchos::RCP<Epetra_Vector> &sat_gas,
                        Teuchos::RCP<Epetra_Vector> &sat_liquid,
                        Teuchos::RCP<Epetra_Vector> &sat_ice) const;
  void DeriveLiquidFlux(const Epetra_Vector &pressure, const Epetra_Vector &rho,
                        const Epetra_Vector &k_rel, const Epetra_Vector &k,
                        const Epetra_Vector &mu,
                        Teuchos::RCP<Epetra_Vector> &flux,
                        double &l2_error) const;
  void DeriveLiquidVelocity(const Epetra_Vector &pressure,
                            const Epetra_Vector &rho,
                            const Epetra_Vector &k_rel,
                            const Epetra_Vector &k,
                            const Epetra_Vector &mu,
                            Teuchos::RCP<Epetra_Vector> &velocity) const;

  // helper function for getting rel perm on the face
  void ComputeUpwindRelPerm(const Epetra_Vector& Pcell,
                            const Epetra_Vector& Pface,
                            const Epetra_Vector& k_rel_cell,
                            const Epetra_Vector& k,
                            const Epetra_Vector& rho,
                            Epetra_Vector& k_rel_face) const;

  // timestepping
  void ComputeUDot(const double t, const Epetra_Vector& u, Epetra_Vector &udot);
  void SetInitialPressureProfileCells(double ref_height, const Epetra_Vector &cell_heights,
                                      Teuchos::RCP<Epetra_Vector> &pressure);
  void SetInitialPressureProfileFaces(double ref_height,
                                      Teuchos::RCP<Epetra_Vector> &pressure);

  // access
  Teuchos::RCP<Epetra_Vector> GetPorosity() { return phi_; }
  Teuchos::RCP<Epetra_Vector> GetTemperature() { return temp_ }
  Teuchos::RCP<Epetra_Vector> GetGasPressure() { return pressure_gas_; }
  Teuchos::RCP<Epetra_Vector> GetGasSaturation() { return sat_gas_; }
  Teuchos::RCP<Epetra_Vector> GetLiquidSaturation() { return sat_liquid_; }
  Teuchos::RCP<Epetra_Vector> GetIceSaturation() { return sat_ice_; }
  Teuchos::RCP<Epetra_Vector> GetPermeability() { return k_; }
  Teuchos::RCP<Epetra_Vector> GetLiquidRelativePermeability() { return k_rel_liquid_; }
  Teuchos::RCP<Epetra_Vector> GetGasDensity() { return rho_gas_; }
  Teuchos::RCP<Epetra_Vector> GetLiquidDensity() { return rho_liquid_; }
  Teuchos::RCP<Epetra_Vector> GetIceDensity() { return rho_ice_; }
  Teuchos::RCP<Epetra_Vector> GetLiquidViscosity() { return mu_liquid_; }
  Teuchos::RCP<Epetra_Vector> GetCellVolumes() { return cell_volumes_; }
  Teuchos::RCP<Epetra_Vector> GetCellHeights() { return cell_heights_; }

  void GetGravity(double g[]) const { for(int i = 0; i < 3; ++i) g[i] = gvec_[i]; }

  Epetra_Vector* CreateCellView(const Epetra_Vector&) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector&) const;

  const Epetra_Map& Map() const { return *dof_map_; }
  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }
  DiffusionMatrix& Matrix() { return *D_; }
  Epetra_Operator& Precon() { return *precon_; }

  // residuals/preconditioners for model evaluator
  //  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time);
  //  void ComputeF(const Epetra_Vector &X, Epetra_Vector &F);

  //  void ComputePrecon(const Epetra_Vector &X);
  //  void ComputePrecon(const Epetra_Vector &X, double h);

private:

  // intermediate steps
  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }
  Teuchos::RCP<Epetra_Map> create_dof_map_(const Epetra_Map&, const Epetra_Map&) const;
  void validate_boundary_conditions_() const;
  DiffusionMatrix* create_diff_matrix_(Teuchos::RCP<AmanziMesh::Mesh>&) const;
  //  void upwind_rel_perm_(const Epetra_Vector&, Epetra_Vector&);

  // maps and access to Vectors
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Epetra_Map> dof_map_;
  Teuchos::RCP<Epetra_Import> face_importer_;
  Teuchos::RCP<Epetra_Import> cell_importer_;
  Teuchos::RCP<Epetra_Vector> cell_volumes_;
  Teuchos::RCP<Epetra_Vector> cell_heights_;

  // constants
  double p_atm_;    // atmospheric pressure
  double gravity_;  // gravitational acceleration (positive coef, directed in -z direction)
  double gvec_[3]; // Vector gravity

  double sigma_gl_; // gas-liquid interfactial tension
  double sigma_il_; // ice-liquid interfactial tension

  const double T0_ = 273.15;
  const double h_iw0;

  bool upwind_k_rel_;

  // constitutive relations
  std::vector< Teuchos::RCP<SaturationCurve> > sat_curves_;
  Teuchos::RCP<EOS> gas_eos_;
  Teuchos::RCP<EOS> liquid_eos_;
  Teuchos::RCP<EOS> ice_eos_;

  // independent variables
  Teuchos::RCP<Epetra_Vector> phi_;
  Teuchos::RCP<Epetra_Vector> temp_;

  // secondary variables
  Teuchos::RCP<Epetra_Vector> pressure_gas_;
  Teuchos::RCP<Epetra_Vector> sat_gas_;
  Teuchos::RCP<Epetra_Vector> sat_liquid_;
  Teuchos::RCP<Epetra_Vector> sat_ice_;

  Teuchos::RCP<Epetra_Vector> rho_gas_;
  Teuchos::RCP<Epetra_Vector> rho_liquid_;
  Teuchos::RCP<Epetra_Vector> rho_ice_;
  Teuchos::RCP<Epetra_Vector> mu_liquid_;

  Teuchos::RCP<Epetra_Vector> k_;
  Teuchos::RCP<Epetra_Vector> k_rel_liquid_;

  // boundary conditions
  Teuchos::RCP<BoundaryFunction> bc_press_;  // Pressure Dirichlet conditions,
                                             //  excluding static head
  Teuchos::RCP<BoundaryFunction> bc_head_;   // Static pressure head conditions,
                                             //  also Dirichlet-type
  Teuchos::RCP<BoundaryFunction> bc_flux_;   // Outward mass flux conditions

  // discretizations, discrete objects
  Teuchos::RCP<MimeticDiscretization> disc_;
  DiffusionPrecon *precon_;
  Teuchos::RCP<DiffusionMatrix> D_;
};

} // close namespace Amanzi

#endif
