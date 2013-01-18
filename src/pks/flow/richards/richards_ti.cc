/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "field_evaluator.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 0

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  ++niter_;

  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);
  double h = t_new - t_old;
  Teuchos::RCP<CompositeVector> u = u_new->data();

  int nc = u->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;
    *out_ << "Richards Residual calculation: T0 = " << t_old << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u)("cell",0,nc) << " " << (*u)("face",0,500) << std::endl;
  }


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
  ApplyDiffusion_(S_next_, res);

  // accumulation term
  AddAccumulation_(res);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    Teuchos::RCP<const CompositeVector> satl1 = S_next_->GetFieldData("saturation_liquid");
    Teuchos::RCP<const CompositeVector> satl0 = S_inter_->GetFieldData("saturation_liquid");
    Teuchos::RCP<const CompositeVector> sati1 = S_next_->GetFieldData("saturation_ice");
    Teuchos::RCP<const CompositeVector> sati0 = S_inter_->GetFieldData("saturation_ice");
    *out_ << "  sat_old_0: " << (*satl0)("cell",0) << ", " << (*sati0)("cell",0) << std::endl;
    *out_ << "  sat_new_0: " << (*satl1)("cell",0) << ", " << (*sati1)("cell",0) << std::endl;
    *out_ << "  sat_old_1: " << (*satl0)("cell",nc) << ", " << (*sati0)("cell",nc) << std::endl;
    *out_ << "  sat_new_1: " << (*satl1)("cell",nc) << ", " << (*sati1)("cell",nc) << std::endl;

    *out_ << "  res0 (after accumulation): " << (*res)("cell",0,0)
          << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after accumulation): " << (*res)("cell",0,nc)
          << " " << (*res)("face",0,500) << std::endl;
  }

#if DEBUG_FLAG
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

  // Dump residual
  int nc = u->data()->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,0) << " "
          << (*u->data())("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u->data())("cell",0,nc) << " "
          << (*u->data())("face",0,500) << std::endl;
  }

  // NaN checking
  Teuchos::RCP<const CompositeVector> pres = u->data();
  for (int c=0; c!=pres->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*pres)("cell",c))) {
      *out_ << "Cutting time step due to NaN in cell residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=pres->size("face",false); ++f) {
    if (boost::math::isnan<double>((*pres)("face",f))) {
      *out_ << "Cutting time step due to NaN in face residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

  // Apply the preconditioner
  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  // Dump correction
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {

  *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " "
        << (*Pu->data())("face",0,3) << std::endl;
  *out_ << "  PC*p1: " << (*Pu->data())("cell",0,nc) << " "
        << (*Pu->data())("face",0,500) << std::endl;
  }

  // NaN checking of correction
  Teuchos::RCP<const CompositeVector> ppres = Pu->data();
  for (int c=0; c!=ppres->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*ppres)("cell",c))) {
      // dump the schur complement
      Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
      std::stringstream filename_s;
      filename_s << "schur_" << S_next_->cycle() << ".txt";
      EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
      *out_ << "updated precon " << S_next_->cycle() << std::endl;

      // print the rel perm
      Teuchos::RCP<const CompositeVector> num_rel_perm =
          S_next_->GetFieldData("numerical_rel_perm");
      Teuchos::RCP<const CompositeVector> rel_perm =
          S_next_->GetFieldData("relative_permeability");
      *out_ << "REL PERM: " << std::endl;
      rel_perm->Print(*out_);
      *out_ << std::endl;
      *out_ << "UPWINDED REL PERM: " << std::endl;
      num_rel_perm->Print(*out_);

      // throw
      *out_ << "Cutting time step due to NaN in PC'd cell residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

  for (int f=0; f!=ppres->size("face",false); ++f) {
    if (boost::math::isnan<double>((*ppres)("face",f))) {
      *out_ << "Cutting time step due to NaN in PC'd face residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Richards::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::RCP<CompositeVector> uw_rel_perm = S_next_->GetFieldData("numerical_rel_perm", name_);

  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "Precon update at t = " << t << std::endl;
  }

  // update state with the solution up.
  S_next_->set_time(t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // Attempt of a hack to deal with zero rel perm
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");
  Teuchos::RCP<CompositeVector> hacked_rel_perm =
      Teuchos::rcp(new CompositeVector(*rel_perm));
  *hacked_rel_perm = *rel_perm;

  double eps = 1.e-12;
  for (int f=0; f!=hacked_rel_perm->size("face"); ++f) {
    if ((*hacked_rel_perm)("face",f) < eps) {
      bc_markers_[f] = Operators::MFD_BC_FLUX;
      bc_values_[f] = 0.0;
      (*hacked_rel_perm)("face",f) = 0.5*eps;
    }
  }

  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");
  
  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->CreateMFDstiffnessMatrices(hacked_rel_perm.ptr());
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_(gvec, hacked_rel_perm, rho, preconditioner_);

  // Update the preconditioner with accumulation terms.
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> dwc_dp =
      S_next_->GetFieldData("dwater_content_d"+key_);
  Teuchos::RCP<const CompositeVector> pres =
      S_next_->GetFieldData(key_);

  // -- update the cell-cell block
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  for (int c=0; c!=dwc_dp->size("cell"); ++c) {
    Acc_cells[c] += (*dwc_dp)("cell",c) / h;
    Fc_cells[c] += (*pres)("cell",c) * (*dwc_dp)("cell",c) / h;

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
  // update the tolerances if we are continuing from an crappy IC
  if (continuation_to_ss_) {
    atol_ = atol0_ + 1.e5*atol0_/(1.0 + S_next_->time());
    rtol_ = rtol0_ + 1.e5*rtol0_/(1.0 + S_next_->time());
  }

  // cell error given by tolerances on water content
  S_next_->GetFieldEvaluator("water_content")->HasFieldChanged(S_next_.ptr(), name_);
  const CompositeVector& wc = *S_next_->GetFieldData("water_content");

  const CompositeVector& res = *du->data();
  double h = S_next_->time() - S_inter_->time();

  double enorm_cell(0.);
  int ncells = res.size("cell");
  for (int c=0; c!=ncells; ++c) {
    double tmp = abs(h*res("cell",c)) / (atol_+rtol_*abs(wc("cell",c)));
    enorm_cell = std::max<double>(enorm_cell, tmp);
  }

  // cell error given by tolerances on pressure
  double enorm_face(0.);
  /*
  int nfaces = res.size("face");
  for (int f=0; f!=nfaces; ++f) {
    double tmp = abs(res("face",f)) / (atol_+rtol_*101325.0);
    enorm_face = std::max<double>(enorm_face, tmp);
  }
  */
  
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    double infnorm_c(0.), infnorm_f(0.);
    res.ViewComponent("cell",false)->NormInf(&infnorm_c);
    res.ViewComponent("face",false)->NormInf(&infnorm_f);

    double buf_c(0.), buf_f(0.);
    MPI_Allreduce(&enorm_cell, &buf_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(&enorm_face, &buf_f, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    Teuchos::OSTab tab = getOSTab();
    *out_ << "ENorm (Infnorm) of: " << name_ << ": "
          << "cell = " << buf_c << " (" << infnorm_c << ")  "
          << "face = " << buf_f << " (" << infnorm_f << ")  " << std::endl;

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



