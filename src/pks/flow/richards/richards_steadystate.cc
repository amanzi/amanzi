#include "EpetraExt_RowMatrixOut.h"
#include "richards_steadystate.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

RegisteredPKFactory<RichardsSteadyState> RichardsSteadyState::reg_("richards steady state");

RichardsSteadyState::RichardsSteadyState(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist,solution),
    Richards(plist, solution) {}

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

  // Create the preconditioner
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");
  mfd_preconditioner_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  mfd_preconditioner_->CreateMFDrhsVectors();

  // update with gravity fluxes
  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), mfd_preconditioner_.ptr());

  // Assemble and precompute the Schur complement for inversion.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (assemble_preconditioner_) {
    mfd_preconditioner_->AssembleGlobalMatrices();
    mfd_preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
    mfd_preconditioner_->UpdatePreconditioner();
  }

  /*
  // dump the schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
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
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  niter_++;

  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->data();

  if (dynamic_mesh_) matrix_->CreateMFDmassMatrices(K_.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << std::setprecision(15);
    *out_ << "----------------------------------------------------------------" << std::endl;
    *out_ << "Residual calculation: t0 = " << t_old
          << " t1 = " << t_new << " h = " << h << std::endl;

    Teuchos::RCP<const CompositeVector> u_old = S_inter_->GetFieldData(key_);

    for (std::vector<int>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      AmanziGeometry::Point c0_centroid = mesh_->cell_centroid(*c0);
      *out_ << "Cell c(" << *c0 << ") centroid = " << c0_centroid << std::endl;

      AmanziMesh::Entity_ID_List fnums0;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

      *out_ << "  p_old(" << *c0 << "): " << (*u_old)("cell",*c0);
      for (int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*u_old)("face",fnums0[n]);
      *out_ << std::endl;

      *out_ << "  p_new(" << *c0 << "): " << (*u)("cell",*c0);
      for (int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*u)("face",fnums0[n]);
      *out_ << std::endl;
    }
  }
#endif

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


    Teuchos::RCP<const CompositeVector> relperm =
        S_next_->GetFieldData("relative_permeability");
    Teuchos::RCP<const Epetra_MultiVector> relperm_bf =
        relperm->ViewComponent("boundary_face",false);
    Teuchos::RCP<const CompositeVector> uw_relperm =
        S_next_->GetFieldData("numerical_rel_perm");
    Teuchos::RCP<const Epetra_MultiVector> uw_relperm_bf =
        uw_relperm->ViewComponent("boundary_face",false);
    Teuchos::RCP<const CompositeVector> flux_dir =
        S_next_->GetFieldData("darcy_flux_direction");

    if (S_next_->HasField("saturation_ice")) {
      Teuchos::RCP<const CompositeVector> sati1 =
          S_next_->GetFieldData("saturation_ice");
      Teuchos::RCP<const CompositeVector> sati0 =
          S_inter_->GetFieldData("saturation_ice");
      for (std::vector<int>::const_iterator c0=dc_.begin();
           c0!=dc_.end(); ++c0) {
        *out_ << "    sat_old(" << *c0 << "): " << (*satl0)("cell",*c0) << ", "
              << (*sati0)("cell",*c0) << std::endl;
        *out_ << "    sat_new(" << *c0 << "): " << (*satl1)("cell",*c0) << ", "
              << (*sati1)("cell",*c0) << std::endl;
      }
    } else {
      for (std::vector<int>::const_iterator c0=dc_.begin();
           c0!=dc_.end(); ++c0) {
        *out_ << "    sat_old(" << *c0 << "): " << (*satl0)("cell",*c0) << std::endl;
        *out_ << "    sat_new(" << *c0 << "): " << (*satl1)("cell",*c0) << std::endl;
      }
    }

    for (std::vector<int>::const_iterator c0=dc_.begin();
         c0!=dc_.end(); ++c0) {
      AmanziMesh::Entity_ID_List fnums0;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

      *out_ << "    k_rel(" << *c0 << "): " << (*uw_relperm)("cell",*c0);
      for (int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*uw_relperm)("face",fnums0[n]);
      *out_ << std::endl;

      *out_ << "  res(" << *c0 << ") (after diffusion): " << (*res)("cell",*c0);
      for (int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*res)("face",fnums0[n]);
      *out_ << std::endl;
    }
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

} // namespace
} // namespace
