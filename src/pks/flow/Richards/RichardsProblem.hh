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
  void Initialize();
  void InitializePressureCells(Teuchos::RCP<Epetra_Vector>& pressure) const;
  void InitializePressureFaces(Teuchos::RCP<Epetra_Vector>& pressure_lambda) const;
  void InitializeSaturation(Teuchos::RCP<Epetra_Vector>& saturation) const;

  // set independent variables
  void SetFluidDensity(double rho);
  void SetFluidViscosity(double mu);
  void SetGravity(const double g[3]);
  void SetGravity(double g); // note this is g in the negative-z-direction

  // derive things
  void DeriveRelPerm(const Epetra_Vector& saturation,
                     Teuchos::RCP<Epetra_Vector>& rel_perm) const;
  void DeriveDarcyFlux(const Epetra_Vector& pressure,
                       const Epetra_Vector& rel_perm,
                       const Epetra_Vector& perm,
                       Teuchos::RCP<Epetra_Vector>& darcy_flux,
                       double *l1_error) const;

  void DeriveDarcyVelocity(const Epetra_Vector& darcy_flux,
                           Teuchos::RCP<Epetra_MultiVector>& darcy_velcoity) const;

private:

  void ValidateBoundaryConditions_() const;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;

  double p_atm_;    // atmospheric pressure
  double rho_;      // fluid density
  double mu_;       // fluid viscosity
  double gravity_;  // gravitational acceleration (positive coef, directed in -z direction)
  Teuchos::RCP<double*> gvec_;  // set to {0,0,-gravity_}

  Teuchos::RCP<Epetra_Vector> cell_volumes_;

  Teuchos::RCP<BoundaryFunction> bc_press_;  // Pressure Dirichlet conditions, excluding static head
  Teuchos::RCP<BoundaryFunction> bc_head_;   // Static pressure head conditions; also Dirichlet-type
  Teuchos::RCP<BoundaryFunction> bc_flux_;   // Outward mass flux conditions

  std::vector< Teuchos::RCP<WaterRetentionBaseModel> > WRM_;

};

} // close namespace Amanzi

#endif
