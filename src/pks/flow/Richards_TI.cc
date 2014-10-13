/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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

#include "LinearOperatorFactory.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void Richards_PK::Functional(double Told, double Tnew, 
                             Teuchos::RCP<CompositeVector> u_old, Teuchos::RCP<CompositeVector> u_new, 
                             Teuchos::RCP<CompositeVector> f)
{ 
  double Tp(Tnew), dTp(Tnew - Told);

  const Epetra_MultiVector& uold_cell = *u_old->ViewComponent("cell");
  const Epetra_MultiVector& unew_cell = *u_new->ViewComponent("cell");

  // update coefficients
  darcy_flux_copy->ScatterMasterToGhosted("face");
  rel_perm_->Compute(*u_new); 
  upwind_->Compute(*darcy_flux_upwind, bc_model, bc_value, *rel_perm_->Krel(), *rel_perm_->Krel(), "k_relative");
  upwind_->Compute(*darcy_flux_upwind, bc_model, bc_value, *rel_perm_->dKdP(), *rel_perm_->dKdP(), "dkdpc");
  UpdateSourceBoundaryData(Tp, *u_new);
  
  // assemble residual for diffusion operator
  op_matrix_->Init();
  op_matrix_->UpdateMatrices(darcy_flux_copy, solution);
  op_matrix_->ApplyBCs();

  Teuchos::RCP<CompositeVector> rhs = op_matrix_->rhs();
  if (src_sink != NULL) AddSourceTerms(*rhs);

  op_matrix_->ComputeNegativeResidual(*u_new, *f);

  // add accumulation term 
  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell");
  Epetra_MultiVector& f_cell = *f->ViewComponent("cell");

  functional_max_norm = 0.0;
  functional_max_cell = 0;

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_->WRM();  
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    double s1, s2, volume;
    for (int i = 0; i < block.size(); i++) {
      int c = block[i];
      s1 = WRM[mb]->saturation(atm_pressure_ - unew_cell[0][c]);
      s2 = WRM[mb]->saturation(atm_pressure_ - uold_cell[0][c]);

      double factor = rho_ * phi[0][c] * mesh_->cell_volume(c) / dTp;
      f_cell[0][c] += (s1 - s2) * factor;

      double tmp = fabs(f_cell[0][c]) / factor;  // calculate errors
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
  op_preconditioner_->ApplyInverse(*X, *Y);
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Richards_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const CompositeVector> u, double dTp)
{
  // update coefficients
  darcy_flux_copy->ScatterMasterToGhosted("face");
  rel_perm_->Compute(*u);

  upwind_->Compute(*darcy_flux_upwind, bc_model, bc_value, *rel_perm_->Krel(), *rel_perm_->Krel(),"k_relative");
  upwind_->Compute(*darcy_flux_upwind, bc_model, bc_value, *rel_perm_->dKdP(), *rel_perm_->dKdP(), "dkdpc");

  // Epetra_MultiVector& dk_face = *rel_perm_->dKdP()->ViewComponent("face");
  // std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_->WRM();
  // const Epetra_IntVector& map_c2mb = rel_perm_->map_c2mb();
  // for (int f = 0; f < nfaces_wghost; f++) {
  //   if (bc_model[f] == Operators::OPERATOR_BC_FACE_NEUMANN && bc_value[f] < 0.0) {
  //     int c = BoundaryFaceGetCell(f);
  //     double face_val = op_matrix_ -> DeriveBoundaryFaceValue(f, *solution, WRM[map_c2mb[c]]);
  //     dk_face[0][f] = -WRM[map_c2mb[c]]->dKdPc (atm_pressure_ - face_val);
  //   }
  //   else if (bc_model[f] == Operators::OPERATOR_BC_FACE_DIRICHLET){
  //     int c = BoundaryFaceGetCell(f);
  //     dk_face[0][f] = -WRM[map_c2mb[c]]->dKdPc (atm_pressure_ - bc_value[f]);
  //   }
  // }

  UpdateSourceBoundaryData(Tp, *u);

  // create diffusion operators
  op_preconditioner_->Init();
  op_preconditioner_->UpdateMatrices(darcy_flux_copy, solution);
  op_preconditioner_->ApplyBCs();

  // add time derivative
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(false);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);

  CompositeVector dSdP(cvs);
  rel_perm_->DerivedSdP(*u->ViewComponent("cell"), *dSdP.ViewComponent("cell"));

  const CompositeVector& phi = *S_->GetFieldData("porosity");
  dSdP.Multiply(rho_, phi, dSdP, 0.0);

  if (dTp > 0.0) {
    op_preconditioner_->AddAccumulationTerm(*u, dSdP, dTp);
  }

  // finalize preconditioner
  int schema_prec_dofs = op_preconditioner_->schema_prec_dofs();
  op_preconditioner_->AssembleMatrix(schema_prec_dofs);
  op_preconditioner_->InitPreconditioner(ti_specs->preconditioner_name, preconditioner_list_); 
}


/* ******************************************************************
* Modify preconditior as needed.
****************************************************************** */
bool Richards_PK::ModifyPredictor(double dT, Teuchos::RCP<const CompositeVector> u0,
                                  Teuchos::RCP<CompositeVector> u)
{
  Teuchos::RCP<CompositeVector> du = Teuchos::rcp(new CompositeVector(*u));
  du->Update(-1.0, *u0, 1.0);
 
  ModifyCorrection(dT, Teuchos::null, u0, du);

  *u = *u0;
  u->Update(1.0, *du, 1.0);
  return false;
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
    const Epetra_IntVector& map = rel_perm_->map_c2mb();

    if (error == buf) {
      int c = functional_max_cell;
      const AmanziGeometry::Point& xp = mesh_->cell_centroid(c);

      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "residual=" << functional_max_norm << " at point";
      for (int i = 0; i < dim; i++) *vo_->os() << " " << xp[i];
      *vo_->os() << std::endl;
 
      c = cell_p;
      const AmanziGeometry::Point& yp = mesh_->cell_centroid(c);

      *vo_->os() << "pressure err=" << error_p << " at point";
      for (int i = 0; i < dim; i++) *vo_->os() << " " << yp[i];
      *vo_->os() << std::endl;

      int mb = map[c];
      double s = (rel_perm_->WRM())[mb]->saturation(atm_pressure_ - uc[0][c]);
      *vo_->os() << "saturation=" << s << " pressure=" << uc[0][c] << std::endl;
    }
  }

  return error;
}


/********************************************************************
* Modifies nonlinear update du based on the maximum allowed change
* of saturation.
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
    Richards_PK::ModifyCorrection(double dT, Teuchos::RCP<const CompositeVector> f,
                                  Teuchos::RCP<const CompositeVector> u,
                                  Teuchos::RCP<CompositeVector> du)
{
  const Epetra_MultiVector& uc = *u->ViewComponent("cell");
  const Epetra_MultiVector& duc = *du->ViewComponent("cell");
  AmanziGeometry::Point face_centr, cell_cntr;
  double max_sat_pert = 0.25;
  double damping_factor = 0.5;
  double reference_pressure = 101325.0;
  
  if (rp_list_.isSublist("clipping parameters")){
    Teuchos::ParameterList& clip_list = rp_list_.sublist("clipping parameters");
    max_sat_pert = clip_list.get<double>("maximum saturation change", 0.25);
    damping_factor = clip_list.get<double>("pressure damping factor", 0.5);
  }

  int nsat_clipped(0), npre_clipped(0);

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_->WRM(); 
  const Epetra_IntVector& map = rel_perm_->map_c2mb();
 
  for (int c = 0; c < ncells_owned; c++) {
    int mb = map[c];
    double pc = atm_pressure_ - uc[0][c];
    double sat = WRM[mb]->saturation(pc);
    double sat_pert;
    if (sat >= 0.5) sat_pert = sat - max_sat_pert;
    else sat_pert = sat + max_sat_pert;
    
    double press_pert = atm_pressure_ - WRM[mb]->capillaryPressure(sat_pert);
    double du_pert_max = fabs(uc[0][c] - press_pert); 
    double tmp = duc[0][c];

    if ((fabs(duc[0][c]) > du_pert_max) && (1 - sat > 1e-5)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "clip saturation: c=" << c 
                   << " p=" << uc[0][c]
                   << " dp: " << duc[0][c] << " -> " << du_pert_max << std::endl;
      }

      if (duc[0][c] >= 0.0) duc[0][c] = du_pert_max;
      else duc[0][c] = -du_pert_max;
      
      nsat_clipped++;
    }    
  }

  for (int c = 0; c < ncells_owned; c++) {
    double unew = uc[0][c] - duc[0][c];
    double tmp = duc[0][c];

    if ((unew > atm_pressure_) && (uc[0][c] < atm_pressure_)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
	 *vo_->os() << "pressure change: " << uc[0][c] << " -> " << unew << std::endl;
      }
      duc[0][c] = tmp * damping_factor;
      npre_clipped++;
    }
  }

  // output statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int nsat_tmp = nsat_clipped, npre_tmp = npre_clipped;
    mesh_->get_comm()->SumAll(&nsat_tmp, &nsat_clipped, 1);
    mesh_->get_comm()->SumAll(&npre_tmp, &npre_clipped, 1);

    if (nsat_clipped > 0 || npre_clipped > 0) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << vo_->color("red") << "saturation/pressure clipped in " 
                 << nsat_clipped << "/" << npre_clipped << " cells" << vo_->reset() << std::endl;
    }
  }

  return (nsat_clipped + npre_clipped) > 0 ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED :
      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

}  // namespace Flow
}  // namespace Amanzi



