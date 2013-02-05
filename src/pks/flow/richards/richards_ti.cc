/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"
#include "Mesh_MSTK.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_RES_FLAG 0

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::fun(double t_old,
                   double t_new,
                   Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new,
                   Teuchos::RCP<TreeVector> g) {
  niter_++;

  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  Teuchos::RCP<CompositeVector> u = u_new->data();

#if DEBUG_FLAG
  int nc = u->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;
    *out_ << "Richards Residual calculation: T0 = " << t_old
          << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3)
          << std::endl;
    *out_ << "  p1: " << (*u)("cell",0,nc) << " " << (*u)("face",0,500)
          << std::endl;
  }
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    Teuchos::RCP<const CompositeVector> satl1 =
        S_next_->GetFieldData("saturation_liquid");
    Teuchos::RCP<const CompositeVector> satl0 =
        S_inter_->GetFieldData("saturation_liquid");

    if (S_next_->HasField("saturation_ice")) {
      Teuchos::RCP<const CompositeVector> sati1 =
          S_next_->GetFieldData("saturation_ice");
      Teuchos::RCP<const CompositeVector> sati0 =
          S_inter_->GetFieldData("saturation_ice");
      *out_ << "    sat_old_0: " << (*satl0)("cell",0) << ", "
            << (*sati0)("cell",0) << std::endl;
      *out_ << "    sat_new_0: " << (*satl1)("cell",0) << ", "
            << (*sati1)("cell",0) << std::endl;
      *out_ << "    sat_old_1: " << (*satl0)("cell",nc) << ", "
            << (*sati0)("cell",nc) << std::endl;
      *out_ << "    sat_new_1: " << (*satl1)("cell",nc) << ", "
            << (*sati1)("cell",nc) << std::endl;
    } else {
      *out_ << "    sat_old_0: " << (*satl0)("cell",0) << std::endl;
      *out_ << "    sat_new_0: " << (*satl1)("cell",0) << std::endl;
      *out_ << "    sat_old_1: " << (*satl0)("cell",nc) << std::endl;
      *out_ << "    sat_new_1: " << (*satl1)("cell",nc) << std::endl;
    }


    *out_ << "  res0 (after diffusion): " << (*res)("cell",0,0)
          << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after diffusion): " << (*res)("cell",0,nc)
          << " " << (*res)("face",0,500) << std::endl;
  }
#endif

  // accumulation term
  AddAccumulation_(res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after accumulation): " << (*res)("cell",0,0)
          << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after accumulation): " << (*res)("cell",0,nc)
          << " " << (*res)("face",0,500) << std::endl;
  }
#endif

#if DEBUG_RES_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << "flow_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << "flow_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif
};

// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void Richards::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

#if DEBUG_FLAG
  // Dump residual
  int nc = u->data()->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,0) << " "
          << (*u->data())("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u->data())("cell",0,nc) << " "
          << (*u->data())("face",0,500) << std::endl;
  }
#endif

  // Apply the preconditioner
  preconditioner_->ApplyInverse(*u->data(), Pu->data().ptr());

#if DEBUG_FLAG
  // Dump correction
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

  *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " "
        << (*Pu->data())("face",0,3) << std::endl;
  *out_ << "  PC*p1: " << (*Pu->data())("cell",0,nc) << " "
        << (*Pu->data())("face",0,500) << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Richards::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon update at t = " << t << std::endl;
  }
#endif

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");
  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");

  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), preconditioner_.ptr());

  // Update the preconditioner with accumulation terms.
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  const Epetra_MultiVector& dwc_dp =
      *S_next_->GetFieldData("dwater_content_d"+key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& pres =
      *S_next_->GetFieldData(key_)->ViewComponent("cell",false);

  // -- update the cell-cell block
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  int ncells = dwc_dp.MyLength();
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += dwc_dp[0][c] / h;
    Fc_cells[c] += pres[0][c] * dwc_dp[0][c] / h;
  }

  // If coupled, update the vector used by surface preconditioner to get the
  // contribution of d overland_source_from_subsurface / d surface_pressure,
  // which is currently in our face unknown.
  //
  // Note this must be done prior to boundary conditions being applied.
  if (coupled_to_surface_via_head_) {
    Epetra_MultiVector& dsource = *S_next_->GetFieldData(
        "doverland_source_from_subsurface_dsurface_pressure", name_)
        ->ViewComponent("cell",false);

    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S_next_->GetMesh("surface"));

    std::vector<Epetra_SerialDenseVector>& Acf_cells = preconditioner_->Acf_cells();

    int ncells_surface = dsource.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- and the cell inside that face
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);

      // -- find my face index in the local numbering
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(cells[0], &faces, &dirs);
      int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();

      // -- set the value
      std::cout << "Setting precon value = " << Acf_cells[cells[0]](my_n) << std::endl;
      dsource[0][c] = Acf_cells[cells[0]](my_n);
    }
  }

  // Assemble and precompute the Schur complement for inversion.
  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (assemble_preconditioner_) {
    preconditioner_->AssembleGlobalMatrices();
    preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
    preconditioner_->UpdatePreconditioner();
  }

  /*
  // dump the schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  *out_ << "updated precon " << S_next_->cycle() << std::endl;

  // print the rel perm
  Teuchos::RCP<const CompositeVector> cell_rel_perm =
      S_next_->GetFieldData("relative_permeability");
  *out_ << "REL PERM: " << std::endl;
  cell_rel_perm->Print(*out_);
  *out_ << std::endl;
  *out_ << "UPWINDED REL PERM: " << std::endl;
  rel_perm->Print(*out_);
  */

};


double Richards::enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  // Calculate water content at the solution.
  S_next_->GetFieldEvaluator("water_content")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wc = *S_next_->GetFieldData("water_content")
      ->ViewComponent("cell",false);

  Teuchos::RCP<const CompositeVector> res = du->data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);
  const Epetra_MultiVector& pres_f = *u->data()->ViewComponent("face",false);
  double h = S_next_->time() - S_inter_->time();

  // Cell error is based upon error in mass conservation relative to
  // the current water content
  double enorm_cell(0.);
  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    double tmp = std::abs(h*res_c[0][c]) / (atol_+rtol_*std::abs(wc[0][c]));
    enorm_cell = std::max<double>(enorm_cell, tmp);
  }

  // Face error given by tolerances on pressure?  Not sure what should be here!
  double enorm_face(0.);
  int nfaces = res_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    double tmp = std::abs(res_f[0][f]) /
        (atol_ + rtol_*( std::abs(pres_f[0][f] - 101325.0) + 101325.));
    enorm_face = std::max<double>(enorm_face, tmp);
  }


  // Write out Inf norms too.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    double infnorm_c(0.), infnorm_f(0.);
    res_c.NormInf(&infnorm_c);
    res_f.NormInf(&infnorm_f);

#ifdef HAVE_MPI
    double buf_c(enorm_cell), buf_f(enorm_face);
    MPI_Allreduce(&buf_c, &enorm_cell, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&buf_f, &enorm_face, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    *out_ << "ENorm (cells) = " << enorm_cell << " (" << infnorm_c << ")  " << std::endl;
    *out_ << "ENorm (faces) = " << enorm_face << " (" << infnorm_f << ")  " << std::endl;
  }

  double enorm_val(std::max<double>(enorm_face, enorm_cell));
#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_val;
};


}  // namespace Flow
}  // namespace Amanzi



