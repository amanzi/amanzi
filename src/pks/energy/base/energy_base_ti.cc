/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Debugger.hh"
#include "boundary_function.hh"
#include "FieldEvaluator.hh"
#include "energy_base.hh"
#include "Op.hh"

namespace Amanzi {
namespace Energy {

#define DEBUG_FLAG 1
#define MORE_DEBUG_FLAG 0

// EnergyBase is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void EnergyBase::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // increment, get timestep
  niter_++;
  double h = t_new - t_old;

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("T_old"); vnames.push_back("T_new");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr()); vecs.push_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);
#endif

  // update boundary conditions
  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  vnames[0] = "K";
  vecs[0] = S_next_->GetFieldData(conductivity_key_).ptr();
  vnames[1] = "uw_K";
  vecs[1] = S_next_->GetFieldData(uw_conductivity_key_).ptr();
  db_->WriteVectors(vnames,vecs,true);

  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  vnames[0] = "e_old";
  vnames[1] = "e_new";
  vecs[0] = S_inter_->GetFieldData(energy_key_).ptr();
  vecs[1] = S_next_->GetFieldData(energy_key_).ptr();
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (acc)", res.ptr());
#endif

  // advection term
  if (implicit_advection_) {
    AddAdvection_(S_next_.ptr(), res.ptr(), true);
  } else {
    AddAdvection_(S_inter_.ptr(), res.ptr(), true);
  }
#if DEBUG_FLAG
  db_->WriteVector("res (adv)", res.ptr());
#endif

  // source terms
  AddSources_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (src)", res.ptr());
#endif

  // Dump residual to state for visual debugging.
#if MORE_DEBUG_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void EnergyBase::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if DEBUG_FLAG
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;
  db_->WriteVector("T_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*T_res", Pu->Data().ptr(), true);
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void EnergyBase::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(*up, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  preconditioner_->Init();
  preconditioner_diff_->Setup(conductivity, Teuchos::null);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  S_next_->GetFieldEvaluator(energy_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& de_dT = *S_next_->GetFieldData(de_dT_key_)
      ->ViewComponent("cell",false);

#if DEBUG_FLAG
  db_->WriteVector("    de_dT", S_next_->GetFieldData(de_dT_key_).ptr());
#endif

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = preconditioner_acc_->local_matrices()->vals;

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

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(S_next_.ptr(), h);

  // update with advection terms
  if (implicit_advection_ && implicit_advection_in_pc_) {
    Teuchos::RCP<const CompositeVector> darcy_flux = S_next_->GetFieldData("darcy_flux");
    S_next_->GetFieldEvaluator(enthalpy_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    Teuchos::RCP<const CompositeVector> dhdT = S_next_->GetFieldData(denthalpy_key_);
    preconditioner_adv_->Setup(*darcy_flux);
    preconditioner_adv_->UpdateMatrices(*darcy_flux, *dhdT);
    ApplyDirichletBCsToEnthalpy_(S_next_.ptr());
    preconditioner_adv_->ApplyBCs(bc_adv_, true);
  }

  // Apply boundary conditions.
  preconditioner_diff_->ApplyBCs(true);
  if (precon_used_) {
    preconditioner_->AssembleMatrix();
    preconditioner_->InitPreconditioner(plist_->sublist("preconditioner"));
  }
};


double EnergyBase::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // Calculate water content at the solution.
  S_next_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& energy = *S_next_->GetFieldData(energy_key_)
      ->ViewComponent("cell",false);

  // Collect additional data.
  Teuchos::RCP<const CompositeVector> res = du->Data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",false);
  const CompositeVector& temp = *u->Data();
  double h = S_next_->time() - S_inter_->time();

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double enorm_cell(-1.);
  int bad_cell = -1;
  unsigned int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    double tmp = std::abs(h*res_c[0][c]) / (atol_ * cv[0][c] + rtol_* std::abs(energy[0][c]));
    if (tmp > enorm_cell) {
      enorm_cell = tmp;
      bad_cell = c;
    }
  }

  bool is_face = res->HasComponent("face");
  double enorm_face(-1.);
  if (is_face) {
    // Face error is mismatch in flux??
    const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);
    int bad_face = -1;
    unsigned int nfaces = res_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
      //    double tmp = flux_tol_ * std::abs(res_f[0][f]) / (atol_+rtol_*273.15);
      double tmp = flux_tol_ * std::abs(h*res_f[0][f])  / (atol_ * cv[0][cells[0]] + rtol_* std::abs(energy[0][cells[0]]));
      if (tmp > enorm_face) {
        enorm_face = tmp;
        bad_face = f;
      }
    }
  }

  // Write out Inf norms too.
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double infnorm_c(0.);
    res_c.NormInf(&infnorm_c);

    ENorm_t err_c;
    ENorm_t l_err_c;
    l_err_c.value = enorm_cell;
    l_err_c.gid = res_c.Map().GID(bad_cell);

    MPI_Allreduce(&l_err_c, &err_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    *vo_->os() << "ENorm (cells) = " << err_c.value << "[" << err_c.gid << "] (" << infnorm_c << ")" << std::endl;

    if (is_face) {
      const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);
      double infnorm_f(0.);
      res_f.NormInf(&infnorm_f);

      ENorm_t err_f;
      ENorm_t l_err_f;
      l_err_f.value = enorm_face;
      l_err_f.gid = res_f.Map().GID(bad_cell);

      MPI_Allreduce(&l_err_f, &err_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      *vo_->os() << "ENorm (faces) = " << err_f.value << "[" << err_f.gid << "] (" << infnorm_f << ")" << std::endl;
    }
  }

  // Communicate and take the max.
  double enorm_val = enorm_cell;
  if (is_face) {
    enorm_val = std::max<double>(enorm_face, enorm_cell);
  }

  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return enorm_val;
};


} // namespace Energy
} // namespace Amanzi
