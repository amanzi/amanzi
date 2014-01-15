/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for the derived MPC for water coupling between surface and subsurface.

   In this method, a Dirichlet BC is used on the subsurface boundary for
   the operator, but the preconditioner is for the flux system with no
   extra unknowns.  On the surface, the TPFA is used, resulting in a
   subsurface-face-only Schur complement that captures all terms.

   ------------------------------------------------------------------------- */
#include "EpetraExt_RowMatrixOut.h"

#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_flux_coupler.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

// -- Setup data.
void MPCSurfaceSubsurfaceFluxCoupler::setup(const Teuchos::Ptr<State>& S) {
  // -- turn off PC assembly: <Parameter name="assemble preconditioner" type="bool" value="false"/>
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");
  plist_->sublist("PKs").sublist(names[0]).set("assemble preconditioner", false);
  plist_->sublist("PKs").sublist(names[1]).set("assemble preconditioner", false);
  plist_->sublist("PKs").sublist(names[0]).set("coupled to surface via flux", true);
  plist_->sublist("PKs").sublist(names[1]).set("coupled to subsurface via flux", true);
  plist_->sublist("PKs").sublist(names[0]).sublist("Diffusion PC").set("coupled to surface", true);
  plist_->sublist("PKs").sublist(names[1]).sublist("Diffusion PC").set("TPFA", true);

  MPCSurfaceSubsurfaceCoupler::setup(S);

  // get the flux key
  flux_key_ = plist_->get<std::string>("flux key");
  S->RequireField(flux_key_, name_)->SetMesh(surf_mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  niter_ = 0;
}

// -- Setup data.
void MPCSurfaceSubsurfaceFluxCoupler::initialize(const Teuchos::Ptr<State>& S) {
  MPCSurfaceSubsurfaceCoupler::initialize(S);
  S->GetField(flux_key_, name_)->set_initialized();
}

// -------------------------------------------------------------
// Special function evaluation process
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  solution_to_state(u_new, S_next_);

  // evaluate the residual of the surface equation
  Teuchos::RCP<TreeVector> surf_u_new = u_new->SubVector(surf_pk_index_);
  Teuchos::RCP<TreeVector> surf_g = g->SubVector(surf_pk_index_);
  surf_pk_->fun(t_old, t_new, Teuchos::null, surf_u_new, surf_g);

  // The residual of the surface equation provides the flux.  This is the mass
  // imbalance on the surface, and is used in the subsurface PK.
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData(flux_key_, name_);
  *source->ViewComponent("cell",false) = *surf_g->Data()->ViewComponent("cell",false);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  Teuchos::RCP<TreeVector> domain_u_new = u_new->SubVector(domain_pk_index_);
  Teuchos::RCP<TreeVector> domain_g = g->SubVector(domain_pk_index_);
  domain_pk_->fun(t_old, t_new, Teuchos::null, domain_u_new, domain_g);

  // Set the surf cell residual to 0.
  surf_g->Data()->ViewComponent("cell",false)->PutScalar(0.);
}


// -------------------------------------------------------------
// Apply preconditioner to u and returns the result in Pu
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  // Apply the combined preconditioner to the subsurface residual
  PreconApply_(u,Pu);

  // Update surface values.
  PreconUpdateSurfaceCells_(Pu);
  PreconUpdateSurfaceFaces_(u,Pu);

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Preconditioner Application:" << std::endl;
  std::vector<std::string> vnames;
  vnames.push_back("  r_sub");
  vnames.push_back("  PC*r_sub");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(u->SubVector(0)->Data().ptr());
  vecs.push_back(Pu->SubVector(0)->Data().ptr());
  domain_db_->WriteVectors(vnames, vecs, true);

  vnames[0] = "  r_surf";
  vnames[1] = "  PC*r_surf";
  vecs[0] = u->SubVector(1)->Data().ptr();
  vecs[1] = Pu->SubVector(1)->Data().ptr();
  surf_db_->WriteVectors(vnames, vecs, true);
}


void MPCSurfaceSubsurfaceFluxCoupler::PreconApply_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Precon applying the subsurf + surf TPFA operator." << std::endl;

  domain_pk_->precon(u->SubVector(domain_pk_index_), Pu->SubVector(domain_pk_index_));
};


bool MPCSurfaceSubsurfaceFluxCoupler::modify_correction(double h,
        Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> du) {
  bool modified = PreconPostprocess_(res, du);
  if (modified) {
    PreconUpdateSurfaceCells_(du);
    PreconUpdateSurfaceFaces_(res,du);
  }
  return modified;
}


// -------------------------------------------------------------
// Post-processing in the preconditioner takes the corrections from the
// surface cell's parent faces and uses them for the surface cell, and then
// calculates a correction for the surface faces.
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::PreconUpdateSurfaceCells_(
    Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Precon updating surface cells." << std::endl;

  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_index_)->Data();

  int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // Correction applies to both the domain face and the surface cell.  Also
  // correct for drift.
  const Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  const Epetra_MultiVector& domain_p_f = *S_next_->GetFieldData("pressure")
      ->ViewComponent("face",false);
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);
  const Epetra_MultiVector& surf_p_c = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);

  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    surf_Pu_c[0][cs] = domain_Pu_f[0][f] - domain_p_f[0][f] + surf_p_c[0][cs];
  }
}


// -------------------------------------------------------------
// Post-processing in the preconditioner takes the corrections from the
// surface cell's parent faces and uses them for the surface cell, and then
// calculates a correction for the surface faces.
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::PreconUpdateSurfaceFaces_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Precon updating surface faces." << std::endl;

  // update delta faces
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_index_)->Data();
  surf_preconditioner_->UpdateConsistentFaceCorrection(*surf_u, surf_Pu.ptr());
}


// updates the preconditioner
void MPCSurfaceSubsurfaceFluxCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;
  niter_++;

  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);

  if (assemble_preconditioner_) {
    mfd_preconditioner_->AssembleGlobalMatrices();
    mfd_preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
            domain_pk_->bc_values());

    mfd_preconditioner_->UpdatePreconditioner();
  }

  /*
  // TEST
  if (S_next_->cycle() == 151) {
  // Dump the Schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);

  std::cout << "CYCLE 176, ITER " << niter_ << "!!!!!!!!" << std::endl;

  changed_solution();
  Teuchos::RCP<TreeVector> up_nc = Teuchos::rcp_const_cast<TreeVector>(up);
  Teuchos::RCP<TreeVector> up2 = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector(*up));
  fun(S_->time(), S_next_->time(), Teuchos::null, up_nc, f1);

  *up2 = *up;
  int c = 4;
  int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
  (*up_nc->SubVector(domain_pk_index_)->Data())("face",f) =
  (*up_nc->SubVector(domain_pk_index_)->Data())("face",f) + .001;
  (*up_nc->SubVector(surf_pk_index_)->Data())("cell",c) =
  (*up_nc->SubVector(surf_pk_index_)->Data())("cell",c) + .001;
  changed_solution();
  fun(S_->time(), S_next_->time(), Teuchos::null, up_nc, f2);

  std::cout << "DFDP: " << std::endl;
  std::cout << "  p0 = " << (*up2->SubVector(domain_pk_index_)->Data())("face",f);
  std::cout << "  sp0 = " << (*up2->SubVector(surf_pk_index_)->Data())("cell",c) << std::endl;
  std::cout << "  p1 = " << (*up_nc->SubVector(domain_pk_index_)->Data())("face",f);
  std::cout << "  sp1 = " << (*up_nc->SubVector(surf_pk_index_)->Data())("cell",c) << std::endl;
  std::cout << "  f0 = " << (*f1->SubVector(domain_pk_index_)->Data())("face",f) << std::endl;
  std::cout << "  f1 = " << (*f2->SubVector(domain_pk_index_)->Data())("face",f) << std::endl;



  double df_dp = ((*f2->SubVector(domain_pk_index_)->Data())("face",f)
  -(*f1->SubVector(domain_pk_index_)->Data())("face",f)) / .001;
  std::cout << "DFDP = " << df_dp << std::endl;
  }
  */
}


bool MPCSurfaceSubsurfaceFluxCoupler::modify_predictor_for_flux_bc_(double h,
        Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Modifying predictor for Flux BCs" << std::endl;

  // To get this to work, we want to first modify the surface predictor,
  // then call the surface residual function to get the surface flux BCs,
  // then call the subsurface modify_predictor.  The subsurface modify
  // routine will alter the surface face lambdas to make them reasonably
  // consistent with the flux BCs.

  // -- call surface's modify_predictor()
  Teuchos::RCP<TreeVector> surf_u = up->SubVector(surf_pk_index_);
  surf_pk_->modify_predictor(h, surf_u);

  // -- call the surface residual function to get the surface flux BCs
  Teuchos::RCP<TreeVector> surf_g = Teuchos::rcp(new TreeVector(*surf_u));
  surf_pk_->fun(S_inter_->time(), S_next_->time(), Teuchos::null,
                surf_u, surf_g);

  // -- set the flux BCs
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData(flux_key_, name_);
  *source->ViewComponent("cell",false) = *surf_g->Data()->ViewComponent("cell",false);

  // -- call the subsurface modify_predictor(), which uses those BCs
  Teuchos::RCP<TreeVector> domain_u = up->SubVector(domain_pk_index_);
  domain_pk_->modify_predictor(h, domain_u);

  // ensure the surface and subsurface match
  modify_predictor_copy_subsurf_to_surf_(h,up);

  return true;
}



// -----------------------------------------------------------------------------
// Modify the nonlinear predictor to avoid issues.
// This is very delicate code to get the order correct...
// -----------------------------------------------------------------------------
bool MPCSurfaceSubsurfaceFluxCoupler::modify_predictor(double h,
        Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Modifying predictor" << std::endl;

  bool changed(false);
  if (modify_predictor_flux_bc_ ||
      (S_next_->cycle() == 0 && modify_predictor_first_flux_bc_)) {
    changed = modify_predictor_for_flux_bc_(h, up);
  } else {
    changed = MPCSurfaceSubsurfaceCoupler::modify_predictor(h, up);
  }
  return changed;
}


} // namespace
