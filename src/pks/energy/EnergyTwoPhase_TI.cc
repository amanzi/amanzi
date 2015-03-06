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
#include "EnergyTwoPhase_PK.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Computes the non-linear functional g = g(t,u,udot)
****************************************************************** */
void EnergyTwoPhase_PK::Functional(
    double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  double h = t_new - t_old;  // get timestep

  // update BCs and conductivity
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());
  UpdateConductivityData(S_.ptr());

  // assemble residual for diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true);

  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *g->Data());

  // add accumulation term
  double dt = t_new - t_old;

  // update the energy at both the old and new times.
  S_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_.ptr(), passwd_);

  const Epetra_MultiVector& e1 = *S_->GetFieldData(energy_key_)->ViewComponent("cell");
  const Epetra_MultiVector& e0 = *S_->GetFieldData(prev_energy_key_)->ViewComponent("cell");
  Epetra_MultiVector& g_c = *g->Data()->ViewComponent("cell");

  int nsize = g_c.MyLength();
  for (int i = 0; i < nsize; ++i) {
    g_c[0][i] += (e1[0][i] - e0[0][i]) / dt;
  }

  // advection term, implicit by default.
  const CompositeVector& flux = *S_->GetFieldData("darcy_flux");
  op_advection_->UpdateMatrices(flux);
};


/* ******************************************************************
* Update the preconditioner at time t and u = up
****************************************************************** */
void EnergyTwoPhase_PK::UpdatePreconditioner(
    double t, Teuchos::RCP<const TreeVector> up, double dt)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "updating preconditioner, T=" << t << std::endl;
  }

  // update BCs and conductivity
  UpdateSourceBoundaryData(t, t + dt, *up->Data());
  UpdateConductivityData(S_.ptr());

  // assemble residual for diffusion operator
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  op_preconditioner_diff_->ApplyBCs(true);

  // update with accumulation terms
  // update the accumulation derivatives, dE/dT
  S_->GetFieldEvaluator(energy_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, "temperature");
  CompositeVector& dEdT = *S_->GetFieldData("denergy_dtemperature", energy_key_);

  if (dt > 0.0) {
    op_acc_->AddAccumulationTerm(*up->Data().ptr(), dEdT, dt, "cell");
  }

  // finalize preconditioner
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);
};


/* ******************************************************************
* TBW
****************************************************************** */
double EnergyTwoPhase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
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

