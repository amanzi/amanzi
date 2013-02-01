#include "richards_steadystate.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 0

RegisteredPKFactory<RichardsSteadyState> RichardsSteadyState::reg_("richards steady state");

void RichardsSteadyState::setup(const Teuchos::Ptr<State>& S) {
  Teuchos::ParameterList bdf_plist = plist_.sublist("time integrator");
  max_iters_ = bdf_plist.get<int>("max iterations", 10);

  Richards::setup(S);
}

// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void RichardsSteadyState::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon update at t = " << t << std::endl;
  }

  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

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


// RichardsSteadyState is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void RichardsSteadyState::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  ++niter_;

  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);
  double h = t_new - t_old;
  Teuchos::RCP<CompositeVector> u = u_new->data();

  int nc = u->size("cell") - 1;
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;
    *out_ << "RichardsSteadyState Residual calculation: T0 = " << t_old << " T1 = " << t_new << " H = " << h << std::endl;
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

    *out_ << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after diffusion): " << (*res)("cell",0,nc) << " " << (*res)("face",0,500) << std::endl;
  }

#if DEBUG_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << "residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << "solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif
};

} // namespace
} // namespace
