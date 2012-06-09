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
  ComputePreconditionerMFD(u, matrix_, Tp, 0.0, false);  // Calculate only stiffness matrix.
  matrix_->ComputeNegativeResidual(u, f);  // compute A*u - g

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const Epetra_Vector& phi = FS->ref_porosity();

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    AmanziMesh::Entity_ID_List block(ncells);
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
  preconditioner_->ApplyInverse(X, Y);
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Richards_PK::update_precon(double Tp, const Epetra_Vector& u, double dTp, int& ierr)
{
  ComputePreconditionerMFD(u, preconditioner_, Tp, dTp, true);
  ierr = 0;

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
     std::printf("Richards PK: updating preconditioner at T(sec)=%10.5e dT(sec)=%9.4e\n", Tp, dTp);
  }
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods. 
****************************************************************** */
double Richards_PK::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error;
  error = ErrorNormSTOMP(u, du);
  // error = ErrorNormRC1(u, du);

  return error;
}


/* ******************************************************************
* Error control a-la STOMP.
****************************************************************** */
double Richards_PK::ErrorNormSTOMP(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error, error_p, error_s;
  if (error_control_ & FLOW_TI_ERROR_CONTROL_PRESSURE) {
    error_p = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      double tmp = fabs(du[c]) / (fabs(u[c] - atm_pressure) + atm_pressure);
      error_p = std::max<double>(error_p, tmp);
    }
  } else {
    error_p = 0.0;
  }

  if (error_control_ & FLOW_TI_ERROR_CONTROL_SATURATION) {
    Epetra_Vector dSdP(mesh_->cell_map(false));
    DerivedSdP(u, dSdP);

    error_s = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      double tmp = dSdP[c] * fabs(du[c]);
      error_s = std::max<double>(error_s, tmp);
    }
  } else {
    error_s = 0.0;
  }
  

#ifdef HAVE_MPI
  double buf = std::max<double>(error_s, error_p);
  du.Comm().MaxAll(&buf, &error, 1);  // find the global maximum
#endif
  return error;
}


/* ******************************************************************
* Error control from RC1
****************************************************************** */
double Richards_PK::ErrorNormRC1(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_norm = 0.0;
  for (int n = 0; n < u.MyLength(); n++) {
    double tmp = fabs(du[n]) / (absolute_tol + relative_tol * fabs(u[n]));
    error_norm = std::max<double>(error_norm, tmp);
  }
 
#ifdef HAVE_MPI
  double buf = error_norm;
  du.Comm().MaxAll(&buf, &error_norm, 1);  // find the global maximum
#endif
  return  error_norm;
}

}  // namespace AmanziFlow
}  // namespace Amanzi



