/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Interface_BDF2.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* .                                                 
****************************************************************** */
Interface_BDF2::Interface_BDF2(Richards_PK* RPK, 
                               Teuchos::ParameterList& rme_list)
{
  RPK_ = RPK;
  rme_list_ = rme_list;
  absolute_tol = rme_list_.get<double>("Absolute error tolerance", 1.0);
  relative_tol = rme_list_.get<double>("Relative error tolerance", 1e-5); 
}


/* ******************************************************************
* .                                                 
****************************************************************** */
void Interface_BDF2::fun(
    const double t, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f)
{
  RPK_->get_matrix()->computeResidual(u, f);  // compute F(u)

  Epetra_Vector* uc = RPK_->get_FS().create_cell_view(u);  
  Epetra_Vector* udotc = RPK_->get_FS().create_cell_view(udot);
  Epetra_Vector* fc = RPK_->get_FS().create_cell_view(f);

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = RPK_->get_mesh();
  Epetra_Vector dS(mesh->cell_map(false));  // compute S'(p)
  RPK_->derivedSdP(*uc, dS);

  Teuchos::RCP<Epetra_Vector> phi = RPK_->get_FS().get_porosity();
  double rho = RPK_->get_rho();

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    double volume = mesh->cell_volume(c);
    double dS_cell = (*phi)[c] * dS[c] * rho * volume;
    (*fc)[c] += dS_cell * (*udotc)[c];
  }
}


/* ******************************************************************
* .                                                 
****************************************************************** */
void Interface_BDF2::precon(const Epetra_Vector& X, Epetra_Vector& Y)
{
  RPK_->get_matrix()->ApplyInverse(X, Y);
}


/* ******************************************************************
* .                                                 
****************************************************************** */
void Interface_BDF2::update_precon(
    const double t, const Epetra_Vector& up, const double dT, int& errc)
{
  RPK_->get_matrix()->update_ML_preconditioner();
  errc = 0;
}


/* ******************************************************************
* .                                                 
****************************************************************** */
double Interface_BDF2::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_norm = 0.0; 
  for (int n=0; n<u.MyLength(); n++) {
    double tmp = abs(du[n]) / (absolute_tol + relative_tol * abs(u[n]));
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

