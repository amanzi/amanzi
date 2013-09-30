/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)

The routine implements interface to BDFx time integrators.  
*/

#include <algorithm>
#include <string>
#include <vector>

#include "Matrix_TPFA.hh"
#include "Richards_PK.hh"
#include "LinearOperatorFactory.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void Richards_PK::fun(
    double Tp, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f, double dTp)
{ 
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    Epetra_Vector* u_cells = FS->CreateCellView(u);
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

    matrix_ -> ApplyBoundaryConditions(bc_model, bc_values);
    AddGravityFluxes_TPFA( Krel_faces, *Grav_term_faces, bc_model, &*matrix_);
    
    Teuchos::RCP<Epetra_Vector> rhs = matrix_->rhs();
    if (src_sink != NULL) AddSourceTerms(src_sink, *rhs);

    matrix_->ComputeNegativeResidual(*u_cells,  f);  
  }
  else {
    // compute inegative residual f = A*u - rhs
    AssembleMatrixMFD(u, Tp);
    matrix_->ComputeNegativeResidual(u, f);     
  }

  const Epetra_Vector& phi = FS->ref_porosity();

  functional_max_norm = 0.0;
  functional_max_cell = 0;

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    double v, s1, s2, volume;
    for (int i = 0; i < block.size(); i++) {
      int c = block[i];
      v = u[c] - udot[c] * dTp;
      s1 = WRM[mb]->saturation(atm_pressure - u[c]);
      s2 = WRM[mb]->saturation(atm_pressure - v);

      double factor = rho_ * phi[c] * mesh_->cell_volume(c) / dTp;
      f[c] += (s1 - s2) * factor;

      double tmp = fabs(f[c]) / factor;  // calculate errors
      if (tmp > functional_max_norm) {
        functional_max_norm = tmp;
        functional_max_cell = c;        
      }
    }
  }
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
void Richards_PK::precon(const Epetra_Vector& X, Epetra_Vector& Y)
{

 if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
  preconditioner_->ApplyInverse(X, Y);
 }
 else {

   Teuchos::ParameterList plist;
   Teuchos::ParameterList& slist = plist.sublist("gmres");
   slist.set<string>("iterative method", "gmres");
   slist.set<double>("error tolerance", 1e-8);
   slist.set<int>("maximum number of iterations", 100);
   Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
   vlist.set("Verbosity Level", "low");


   AmanziSolvers::LinearOperatorFactory<Matrix_MFD, Epetra_Vector, Epetra_BlockMap> factory;
   Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, Epetra_Vector, Epetra_BlockMap> > 
     solver = factory.Create("gmres", plist, preconditioner_, preconditioner_);
   
 
   solver->ApplyInverse(X,Y);

   
 }
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Richards_PK::update_precon(double Tp, const Epetra_Vector& u, double dTp, int& ierr)
{
  AssemblePreconditionerMFD(u, Tp, dTp);
}


/* ******************************************************************
* .                                   
****************************************************************** */
int Richards_PK::ApllyPrecInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  return preconditioner_->ApplyInverse(X, Y);
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
  double error, error_p, error_r;
  int cell_p, cell_r;

  if (error_control_ & FLOW_TI_ERROR_CONTROL_PRESSURE) {
    error_p = 0.0;
    cell_p = 0;
    for (int c = 0; c < ncells_owned; c++) {
      double tmp = fabs(du[c]) / (fabs(u[c] - atm_pressure) + atm_pressure);
      if (tmp > error_p) {
        error_p = tmp;
        cell_p = c;
      } 
    }
  } else {
    error_p = 0.0;
  }

  if (error_control_ & FLOW_TI_ERROR_CONTROL_RESIDUAL) {
    error_r = functional_max_norm;
  } else {
    error_r = 0.0;
  }

  error = std::max(error_r, error_p);

#ifdef HAVE_MPI
  double buf = error;
  du.Comm().MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  // maximum error is printed out only on one processor
  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    if (error == buf) {
      int c = functional_max_cell;
      const AmanziGeometry::Point& xp = mesh_->cell_centroid(c);

      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "residual=" << functional_max_norm << " at point";
      for (int i = 0; i < dim; i++) *(vo_->os()) << " " << xp[i];
      *(vo_->os()) << endl;
 
      c = cell_p;
      const AmanziGeometry::Point& yp = mesh_->cell_centroid(c);

      *(vo_->os()) << "pressure err=" << error_p << " at point";
      for (int i = 0; i < dim; i++) *(vo_->os()) << " " << yp[i];
      *(vo_->os()) << endl;

      int mb = (rel_perm->map_c2mb())[c];
      double s = (rel_perm->WRM())[mb]->saturation(atm_pressure - u[c]);
      *(vo_->os()) << "saturation=" << s << " pressure=" << u[c] << endl;
    }
  }

  // if (error_control_ & FLOW_TI_ERROR_CONTROL_SATURATION) {
  // }
  
  return error;
}


/* ******************************************************************
* Error control from RC1 (OBSOLETE)
****************************************************************** */
double Richards_PK::ErrorNormRC1(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_norm = 0.0;
  double absolute_tol = 1.0, relative_tol = 1e-6;
  for (int n = 0; n < u.MyLength(); n++) {
    double tmp = fabs(du[n]) / (absolute_tol + relative_tol * fabs(u[n]));
    error_norm = std::max(error_norm, tmp);
  }
 
#ifdef HAVE_MPI
  double buf = error_norm;
  du.Comm().MaxAll(&buf, &error_norm, 1);  // find the global maximum
#endif
  return  error_norm;
}


/********************************************************************
* Modifies nonlinear update du based on the maximum allowed change
* of saturation.
****************************************************************** */
bool Richards_PK::modify_update_step(double h, Epetra_Vector& u, Epetra_Vector& du)
{
  double max_sat_pert = 0.25;
  bool ret_val = false;
  double dumping_factor = 0.6;
  double reference_pressure = 101325.;

  int ncells_clipped(0);
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM(); 
 
  for (int c = 0; c < ncells_owned; c++) {
    int mb = (rel_perm->map_c2mb())[c];
    double pc = atm_pressure - u[c];
    double sat = WRM[mb]->saturation(pc);
    double sat_pert;
    if (sat >= 0.5) sat_pert = sat - max_sat_pert;
    else sat_pert = sat + max_sat_pert;
    
    double press_pert = atm_pressure - WRM[mb]->capillaryPressure(sat_pert);
    double du_pert_max = fabs(u[c] - press_pert); 

    if (fabs(du[c]) > du_pert_max) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << "saturation clipping in cell " << c 
                     << " pressure change: " << du[c] << " -> " << du_pert_max << endl;
      }

      double tmp = du[c];
       
      if (du[c] >= 0.0) du[c] = fabs(du_pert_max);
      else du[c] = -fabs(du_pert_max);

      ncells_clipped++;
      ret_val = true;
    }    
  }

  //  if (verbosity >= FLOW_VERBOSITY_HIGH) {
  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    int ncells_tmp = ncells_clipped;
    du.Comm().SumAll(&ncells_tmp, &ncells_clipped, 1);
    if (MyPID == 0 && ncells_clipped > 0) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "saturation was clipped in " << ncells_clipped << " cells" << endl;
    }
  }

  return ret_val;
}


/********************************************************************
* Converts the BDF1 time intgerator to the Newton solver
****************************************************************** */
bool Richards_PK::IsPureNewton() const
{
  if (experimental_solver_ == FLOW_SOLVER_NKA) return false;
  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) return false;
  return true; 
}

}  // namespace AmanziFlow
}  // namespace Amanzi



