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

#define DEBUG_FLAG 1

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);
  double h = t_new - t_old;
  Teuchos::RCP<CompositeVector> u = u_new->data();

  int nc = u->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;
    *out_ << "Richards Residual calculation: T0 = " << t_old << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u)("cell",0,nc) << " " << (*u)("face",0,497) << std::endl;
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

    *out_ << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after diffusion): " << (*res)("cell",0,nc) << " " << (*res)("face",0,497) << std::endl;
  }

  // accumulation term
  AddAccumulation_(res);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after accumulation): " << (*res)("cell",0,0)
          << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after accumulation): " << (*res)("cell",0,nc)
          << " " << (*res)("face",0,497) << std::endl;
  }
};

// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void Richards::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  // Dump residual
  int nc = u->data()->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,0) << " "
          << (*u->data())("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u->data())("cell",0,nc) << " "
          << (*u->data())("face",0,497) << std::endl;
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

  *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " "
        << (*Pu->data())("face",0,3) << std::endl;
  *out_ << "  PC*p1: " << (*Pu->data())("cell",0,nc) << " "
        << (*Pu->data())("face",0,497) << std::endl;
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
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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
  UpdateBoundaryConditionsPreconditioner_();

  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");
  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");

  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_(gvec, rel_perm, rho, preconditioner_);

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

}  // namespace Flow
}  // namespace Amanzi



