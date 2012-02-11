/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)

The routine implement interface to BDFx time integrators.  
*/

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate f(u, du/dt) = s du/dt + A*u - g.                                         
****************************************************************** */
void Richards_PK::fun(double T, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f)
{
  T_internal = T_physical = T;
  computePreconditionerMFD(u, matrix, T, 0.0, false);  // Calculate only stiffness matrix.
  matrix->computeNegativeResidual(u, f);  // compute A*u - g

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const Epetra_Vector& phi = FS->ref_porosity();

  Epetra_Vector* u_cells = FS->createCellView(u);
  Epetra_Vector dSdP(mesh_->cell_map(false));
  derivedSdP(*u_cells, dSdP);

  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho * phi[c] * dSdP[c] * volume;
    f[c] += factor * udot[c]; 
  }
}


/* ******************************************************************
* .                                                 
****************************************************************** */
void Richards_PK::precon(const Epetra_Vector& X, Epetra_Vector& Y)
{
  preconditioner->ApplyInverse(X, Y);
}


/* ******************************************************************
* Compute new preconditioner B(p, dT_prec). For BDF2 method, we need
* a separate memory allocation.                                              
****************************************************************** */
void Richards_PK::update_precon(double Tp, const Epetra_Vector& u, double dTp, int& ierr)
{
  computePreconditionerMFD(u, preconditioner, Tp, dTp, true);
  ierr = 0;
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.                                                 
****************************************************************** */
double Richards_PK::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_norm = 0.0; 
  for (int n=0; n<u.MyLength(); n++) {
    double tmp = abs(du[n]) / (absolute_tol_sss + relative_tol_sss * abs(u[n]));
    error_norm = std::max<double>(error_norm, tmp);
  }

  // find the global maximum
#ifdef HAVE_MPI
  double buf = error_norm;
  //MPI_Allreduce(&buf, &error_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  du.Comm().MaxAll(&buf, &error_norm, 1);
#endif
  return  error_norm;
}


}  // namespace AmanziFlow
}  // namespace Amanzi



