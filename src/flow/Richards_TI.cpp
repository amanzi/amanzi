/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* .                                                 
****************************************************************** */
void Richards_PK::fun(
    const double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f)
{
  matrix->computeResidual(u, f);  // compute F(u)

  Epetra_Vector* u_cells = FS->createCellView(u);  
  Epetra_Vector dS(mesh_->cell_map(false));
  derivedSdP(*u_cells, dS);

  Teuchos::RCP<Epetra_Vector> phi = FS->get_porosity();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double dS_cell = (*phi)[c] * dS[c] * rho * volume;
    f[c] -= dS_cell * u[c];
  }
}


/* ******************************************************************
* .                                                 
****************************************************************** */
void Richards_PK::precon(const Epetra_Vector& X, Epetra_Vector& Y)
{
  matrix->ApplyInverse(X, Y);
}


/* ******************************************************************
* Compute new preconsitioner Sff(p, dT).                                                 
****************************************************************** */
void Richards_PK::update_precon(
    const double T, const Epetra_Vector& u, const double dT_bdf2, int& ierr)
{
  dT = dT_bdf2;  // this dT will be used in calculation of preconditioner 
  computePreconditionerMFD(u, matrix);
  ierr = 0;
}


/* ******************************************************************
* .                                                 
****************************************************************** */
double Richards_PK::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_norm = 0.0; 
  for (int n=0; n<u.MyLength(); n++) {
    double tmp = abs(du[n]) / (absolute_tol_bdf + relative_tol_bdf * abs(u[n]));
    error_norm = std::max<double>(error_norm, tmp);
  }

  // find the global maximum
#ifdef HAVE_MPI
  double buf = error_norm;
  MPI_Allreduce(&buf, &error_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return  error_norm;
}

}  // namespace AmanziFlow
}  // namespace Amanzi



