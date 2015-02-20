/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  EnergyBase is a BDFFnBase
*/

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "boundary_function.hh"
#include "FieldEvaluator.hh"
#include "Energy_PK.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Computes the non-linear functional g = g(t,u,udot)
****************************************************************** */
void Energy_PK::Functional(
    double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  double h = t_new - t_old;  // get timestep

  // pointer-copy temperature into states and update any auxilary data
  // solution_to_state(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  UpdateSourceBoundaryData(t_old, t_new, *u);

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  // ApplyDiffusion_(S_next_.ptr(), res.ptr());

  // accumulation term
  // AddAccumulation_(res.ptr());

  // advection term, implicit by default, options for explicit
  // AddAdvection_(S_next_.ptr(), res.ptr(), true);

  // source terms
  // AddSources_(S_next_.ptr(), res.ptr());
};


/* ******************************************************************
* Apply the preconditioner to u and return the result in Pu.
****************************************************************** */
void Energy_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                    Teuchos::RCP<TreeVector> Pu)
{
  op_preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
}


/* ******************************************************************
* Update the preconditioner at time t and u = up
****************************************************************** */
void Energy_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Precon update at t = " << t << std::endl;
  }

  // update state with the solution up.
  // PKDefaultBase::solution_to_state(*up, S_next_);

  // update boundary conditions
  bc_temperature->Compute(S_->time());
  bc_flux->Compute(S_->time());
  UpdateSourceBoundaryData(t, t + h, *up->Data());

  // div K_e grad u
  UpdateConductivityData(S_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity = S_->GetFieldData(uw_conductivity_key_);

  // assemble residual for diffusion operator
  op_matrix_->Init();
  // op_matrix_diff_->Setup(conductivity.ptr());
  // op_matrix_diff_->UpdateMatrices(darcy_flux_copy, *up->Data());
  op_matrix_diff_->ApplyBCs(op_bc_);
  // op_matrix_->ComputeNegativeResidual(*u_new, *f);

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  /*
  S_->GetFieldEvaluator(energy_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, key_);
  const Epetra_MultiVector& de_dT = *S_next_->GetFieldData(de_dT_key_)->ViewComponent("cell",false);

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();

  // -- update the diagonal
  unsigned int ncells = de_dT.MyLength();

  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {
    // do not add in de/dT if the height is 0
    const Epetra_MultiVector& pres = *S_next_->GetFieldData("surface_pressure")
        ->ViewComponent("cell",false);
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
    for (unsigned int c=0; c!=ncells; ++c) {
      Acc_cells[c] += pres[0][c] >= patm ? de_dT[0][c] / h : 0.;
    }
  } else {
    for (unsigned int c=0; c!=ncells; ++c) {
      Acc_cells[c] += de_dT[0][c] / h;
    }
  }
  */

  // -- update preconditioner with source term derivatives if needed
  // AddSourcesToPrecon_(S_next_.ptr(), h);

  // Apply boundary conditions.
  op_preconditioner_diff_->ApplyBCs(op_bc_);
};


/* ******************************************************************
* TBW
****************************************************************** */
double Energy_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                            Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Calculate water content at the solution.
  S_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& energy = *S_->GetFieldData(energy_key_)->ViewComponent("cell", false);

  // Collect additional data.
  Teuchos::RCP<const CompositeVector> res = du->Data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);
  /*
  const Epetra_MultiVector& cv = *S_->GetFieldData(cell_vol_key_)->ViewComponent("cell", false);
  const CompositeVector& temp = *u->Data();
  double h = S_->time() - S_->time();

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double enorm_cell(-1.);
  int bad_cell = -1;
  unsigned int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    double tmp = std::abs(h*res_c[0][c]) / (atol_ * cv[0][c]*2.e6 + rtol_* std::abs(energy[0][c]));
    if (tmp > enorm_cell) {
      enorm_cell = tmp;
      bad_cell = c;
    }
  }

  // Face error is mismatch in flux??
  double enorm_face(-1.);
  int bad_face = -1;
  unsigned int nfaces = res_f.MyLength();
  for (unsigned int f=0; f!=nfaces; ++f) {
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
    //    double tmp = flux_tol_ * std::abs(res_f[0][f]) / (atol_+rtol_*273.15);
    double tmp = flux_tol_ * std::abs(h*res_f[0][f]) / 
        (atol_ * cv[0][cells[0]]*2.e6 + rtol_* std::abs(energy[0][cells[0]]));
    if (tmp > enorm_face) {
      enorm_face = tmp;
      bad_face = f;
    }
  }

  // Write out Inf norms too.
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double infnorm_c(0.), infnorm_f(0.);
    res_c.NormInf(&infnorm_c);
    res_f.NormInf(&infnorm_f);

    ENorm_t err_f, err_c;
#ifdef HAVE_MPI
    ENorm_t l_err_f, l_err_c;
    l_err_f.value = enorm_face;
    l_err_f.gid = res_f.Map().GID(bad_face);
    l_err_c.value = enorm_cell;
    l_err_c.gid = res_c.Map().GID(bad_cell);

    MPI_Allreduce(&l_err_c, &err_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    MPI_Allreduce(&l_err_f, &err_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#else
    err_f.value = enorm_face;
    err_f.gid = bad_face;
    err_c.value = enorm_cell;
    err_c.gid = bad_cell;
#endif

    *vo_->os() << "ENorm (cells) = " << err_c.value << "[" << err_c.gid << "] (" << infnorm_c << ")" << std::endl;
    *vo_->os() << "ENorm (faces) = " << err_f.value << "[" << err_f.gid << "] (" << infnorm_f << ")" << std::endl;
  }

  // Communicate and take the max.
  double enorm_val(std::max<double>(enorm_face, enorm_cell));
#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_val;
  */
};

}  // namespace Energy
}  // namespace Amanzi

