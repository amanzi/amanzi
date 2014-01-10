/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The routine implements interface to the BDF1 time integrator.  
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
void Richards_PK::Functional(double Told, double Tnew, 
                             Teuchos::RCP<CompositeVector> u_old, Teuchos::RCP<CompositeVector> u_new, 
                             Teuchos::RCP<CompositeVector> f)
{ 
  double Tp(Tnew), dTp(Tnew - Told);

  const Epetra_MultiVector& uold_cells = *u_old->ViewComponent("cell");
  const Epetra_MultiVector& unew_cells = *u_new->ViewComponent("cell");

  AssembleMatrixMFD(*u_new, Tp);
  matrix_->ComputeNegativeResidual(*u_new, *f);

  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell");
  Epetra_MultiVector& f_cells = *f->ViewComponent("cell");



  functional_max_norm = 0.0;
  functional_max_cell = 0;

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    double s1, s2, volume;
    for (int i = 0; i < block.size(); i++) {
      int c = block[i];
      s1 = WRM[mb]->saturation(atm_pressure_ - unew_cells[0][c]);
      s2 = WRM[mb]->saturation(atm_pressure_ - uold_cells[0][c]);

      double factor = rho_ * phi[0][c] * mesh_->cell_volume(c) / dTp;
      f_cells[0][c] += (s1 - s2) * factor;

      double tmp = fabs(f_cells[0][c]) / factor;  // calculate errors
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
void Richards_PK::ApplyPreconditioner(Teuchos::RCP<const CompositeVector> X, 
                                      Teuchos::RCP<CompositeVector> Y)
{
  preconditioner_->ApplyPreconditioner(*X, *Y);
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Richards_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const CompositeVector> u, double dTp)
{
  AssemblePreconditionerMFD(*u, Tp, dTp);
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods. 
****************************************************************** */
double Richards_PK::ErrorNorm(Teuchos::RCP<const CompositeVector> u, 
                              Teuchos::RCP<const CompositeVector> du)
{
  double error;
  error = ErrorNormSTOMP(*u, *du);

  return error;
}


/* ******************************************************************
* Error control a-la STOMP.
****************************************************************** */
double Richards_PK::ErrorNormSTOMP(const CompositeVector& u, const CompositeVector& du)
{
  const Epetra_MultiVector& uc = *u.ViewComponent("cell");
  const Epetra_MultiVector& duc = *du.ViewComponent("cell");

  double error, error_p, error_r;
  int cell_p(0), cell_r(0);

  if (error_control_ & FLOW_TI_ERROR_CONTROL_PRESSURE) {
    error_p = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - atm_pressure_) + atm_pressure_);
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
    const Epetra_MultiVector& map_c2mb = *rel_perm->map_c2mb().ViewComponent("cell");

    if (error == buf) {
      int c = functional_max_cell;
      const AmanziGeometry::Point& xp = mesh_->cell_centroid(c);

      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "residual=" << functional_max_norm << " at point";
      for (int i = 0; i < dim; i++) *vo_->os() << " " << xp[i];
      *vo_->os() << endl;
 
      c = cell_p;
      const AmanziGeometry::Point& yp = mesh_->cell_centroid(c);

      *vo_->os() << "pressure err=" << error_p << " at point";
      for (int i = 0; i < dim; i++) *vo_->os() << " " << yp[i];
      *vo_->os() << endl;

      int mb = map_c2mb[0][c];
      double s = (rel_perm->WRM())[mb]->saturation(atm_pressure_ - uc[0][c]);
      *vo_->os() << "saturation=" << s << " pressure=" << uc[0][c] << endl;
    }
  }

  return error;
}


/********************************************************************
* Modifies nonlinear update du based on the maximum allowed change
* of saturation.
****************************************************************** */
bool Richards_PK::ModifyCorrection(
    double dT, Teuchos::RCP<const CompositeVector> f,
    Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> du)
{
  const Epetra_MultiVector& uc = *u->ViewComponent("cell");
  const Epetra_MultiVector& duc = *du->ViewComponent("cell");

  double max_sat_pert = 0.25;
  double damping_factor = 0.6;
  double reference_pressure = 101325.0;

  int ncells_clipped(0);

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM(); 
  const Epetra_MultiVector& map_c2mb = *rel_perm->map_c2mb().ViewComponent("cell");
 
  for (int c = 0; c < ncells_owned; c++) {
    int mb = map_c2mb[0][c];
    double pc = atm_pressure_ - uc[0][c];
    double sat = WRM[mb]->saturation(pc);
    double sat_pert;
    if (sat >= 0.5) sat_pert = sat - max_sat_pert;
    else sat_pert = sat + max_sat_pert;
    
    double press_pert = atm_pressure_ - WRM[mb]->capillaryPressure(sat_pert);
    double du_pert_max = fabs(uc[0][c] - press_pert); 

    if ((fabs(duc[0][c]) > du_pert_max) && (1 - sat > 1e-5)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "clip saturation: c=" << c 
                   << " p=" << uc[0][c]
                   << " dp: " << duc[0][c] << " -> " << du_pert_max << endl;
      }

      if (duc[0][c] >= 0.0) duc[0][c] = du_pert_max;
      else duc[0][c] = -du_pert_max;
      
      ncells_clipped++;
    }    
  }

  for (int c = 0; c < ncells_owned; c++) {
    double unew = uc[0][c] - duc[0][c];
    double tmp = duc[0][c];

    if ((unew > atm_pressure_) && (uc[0][c] < atm_pressure_)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
	 *vo_->os() << "pressure change: " << uc[0][c] << " -> " << unew << endl;
      }
      duc[0][c] = tmp * damping_factor;
      ncells_clipped++;
    }
  }

  // output statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int ncells_tmp = ncells_clipped;
    mesh_->get_comm()->SumAll(&ncells_tmp, &ncells_clipped, 1);

    if (MyPID == 0 && ncells_clipped > 0) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << vo_->color("red") << "saturation was clipped in " 
                 << ncells_clipped << " cells" << vo_->reset() << endl;
    }
  }

  return (ncells_clipped > 0);
}

}  // namespace AmanziFlow
}  // namespace Amanzi



