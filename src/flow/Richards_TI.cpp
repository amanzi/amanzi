/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)

The routine implement interface to BDFx time integrators.  
*/

#include <algorithm>
#include <string>
#include <vector>

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate f(u, du/dt) = d(s(u))/dt + A*u - g.                                         
****************************************************************** */
void Richards_PK::fun(
    double Tp, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f, double dTp)
{
  // T_internal = Tp;  breaks internal clock (lipnikov@lanl.gov)
  ComputePreconditionerMFD(u, matrix, Tp, 0.0, false);  // Calculate only stiffness matrix.
  matrix->computeNegativeResidual(u, f);  // compute A*u - g

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const Epetra_Vector& phi = FS->ref_porosity();

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    double v, s1, s2, volume;
    for (int i = 0; i < block.size(); i++) {
      int c = block[i];
      v = u[c] - udot[c] * dTp;
      s1 = WRM[mb]->saturation(atm_pressure - u[c]);
      s2 = WRM[mb]->saturation(atm_pressure - v);

      volume = mesh_->cell_volume(c);
      f[c] += rho * phi[c] * volume * (s1 - s2) / dTp;
    }
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
  ComputePreconditionerMFD(u, preconditioner, Tp, dTp, true);
  ierr = 0;
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.                                                 
****************************************************************** */
double Richards_PK::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_norm = 0.0;
  for (int n = 0; n < u.MyLength(); n++) {
    double tmp = fabs(du[n]) / (absolute_tol + relative_tol * fabs(u[n]));
    error_norm = std::max<double>(error_norm, tmp);
  }

  // find the global maximum
#ifdef HAVE_MPI
  double buf = error_norm;
  // MPI_Allreduce(&buf, &error_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  du.Comm().MaxAll(&buf, &error_norm, 1);
#endif
  return  error_norm;
}


}  // namespace AmanziFlow
}  // namespace Amanzi



