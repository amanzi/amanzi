/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

MPC for the Coupled Permafrost model.  This MPC sits at the top of the
subtree:

                    MPCPermafrost
                     /          \
                    /            \
                   /              \
         surf/subsurf            surf/subsurf
           water                   energy
         /      \                  /      \
        /        \                /        \
    flow/        flow/         energy/     energy/
  permafrost  icy_overland    threephase    surface_ice

------------------------------------------------------------------------- */
#include "FieldEvaluator.hh"

#include "permafrost_model.hh"
#include "MatrixMFD_Coupled_Surf.hh"
#include "MatrixMFD_Surf.hh"
#include "mpc_delegate_ewc.hh"
#include "mpc_surface_subsurface_flux_coupler.hh"
#include "mpc_permafrost.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void MPCPermafrost::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::setup(S);

  // ensure the order of PKs is correct
  // TODO --etc

  // option to decoupled and revert to StrongMPC
  decoupled_ = plist_->get<bool>("decoupled",false);

  // off-diagonal terms needed by MPCCoupledCells
  A_key_ = plist_->get<std::string>("conserved quantity A", "water_content");
  B_key_ = plist_->get<std::string>("conserved quantity B", "energy");
  y1_key_ = plist_->get<std::string>("primary variable A", "pressure");
  y2_key_ = plist_->get<std::string>("primary variable B", "temperature");
  dA_dy2_key_ = std::string("d")+A_key_+std::string("_d")+y2_key_;
  dB_dy1_key_ = std::string("d")+B_key_+std::string("_d")+y1_key_;

  Key mesh_key = plist_->get<std::string>("mesh key", "domain");
  domain_mesh_ = S->GetMesh(mesh_key);

  Key surf_mesh_key = plist_->get<std::string>("surface mesh key", "surface");
  surf_mesh_ = S->GetMesh(surf_mesh_key);

  // preconditioner
  Teuchos::ParameterList pc_sublist = plist_->sublist("Coupled PC");
  mfd_surf_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Coupled_Surf(
      pc_sublist, domain_mesh_));

  // Set the subblocks.  Note these are the flux-coupled PCs, which are
  // MatrixMFD_Surfs.
  Teuchos::RCP<Operators::MatrixMFD_Surf> pcA = sub_pks_[0]->coupled_preconditioner();
  Teuchos::RCP<Operators::MatrixMFD_Surf> pcB = sub_pks_[1]->coupled_preconditioner();

  Teuchos::RCP<Operators::MatrixMFD_TPFA> pcA_surf = pcA->GetSurfaceOperator();
  Teuchos::RCP<Operators::MatrixMFD_TPFA> pcB_surf = pcB->GetSurfaceOperator();

  mfd_surf_preconditioner_->SetSubBlocks(pcA, pcB);
  mfd_surf_preconditioner_->SetSurfaceOperators(pcA_surf, pcB_surf);

  // setup and initialize the preconditioner
  mfd_surf_preconditioner_->SymbolicAssembleGlobalMatrices();
  mfd_surf_preconditioner_->InitPreconditioner();

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_->get<std::string>("predictor type", "none");
  if (predictor_string == "none") {
    predictor_type_ = PREDICTOR_NONE;
  } else if (predictor_string == "ewc") {
    predictor_type_ = PREDICTOR_EWC;
  } else if (predictor_string == "smart ewc") {
    predictor_type_ = PREDICTOR_SMART_EWC;
  } else {
    Errors::Message message(std::string("Invalid predictor type ")+predictor_string);
    Exceptions::amanzi_throw(message);
  }

  // create the EWC delegate if requested.
  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    ewc_ = Teuchos::rcp(new MPCDelegateEWC(*plist_));

    Teuchos::RCP<PermafrostModel> model = Teuchos::rcp(new PermafrostModel());
    ewc_->set_model(model);
    ewc_->setup(S);
  }

  // grab the PKs
  coupled_flow_pk_ = sub_pks_[0];
  coupled_energy_pk_ = sub_pks_[1];

  // grab the debuggers
  domain_db_ = coupled_flow_pk_->domain_debugger();
  surf_db_ = coupled_flow_pk_->surface_debugger();
}


void MPCPermafrost::initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::initialize(S);
  if (ewc_ != Teuchos::null) ewc_->initialize(S);
}


void MPCPermafrost::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::set_states(S,S_inter,S_next);
  if (ewc_ != Teuchos::null) ewc_->set_states(S,S_inter,S_next);
}


void MPCPermafrost::commit_state(double dt, const Teuchos::RCP<State>& S) {
  StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::commit_state(dt,S);
  if (ewc_ != Teuchos::null) ewc_->commit_state(dt,S);
}


// update the predictor to be physically consistent
bool MPCPermafrost::modify_predictor(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Modifying predictor, MPCPermafrost." << std::endl;

    std::vector<std::string> vnames(2);
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs(2);

    vnames[0] = "surf p";
    vnames[1] = "surf T";
    vecs[0] = up->SubVector(0)->SubVector(1)->Data().ptr();
    vecs[1] = up->SubVector(1)->SubVector(1)->Data().ptr();
    surf_db_->WriteVectors(vnames, vecs, true);

    vnames[0] = "sub p";
    vnames[1] = "sub T";
    vecs[0] = up->SubVector(0)->SubVector(0)->Data().ptr();
    vecs[1] = up->SubVector(1)->SubVector(0)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  bool changed(false);
  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    // make a new TreeVector that is just the subsurface values (by pointer).
    Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector());
    domain_u_tv->PushBack(up->SubVector(0)->SubVector(0));
    domain_u_tv->PushBack(up->SubVector(1)->SubVector(0));
    changed = ewc_->modify_predictor(h, domain_u_tv);

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Modifying predictor, EWC." << std::endl;

      std::vector<std::string> vnames(2);
      std::vector< Teuchos::Ptr<const CompositeVector> > vecs(2);

      vnames[0] = "surf p";
      vnames[1] = "surf T";
      vecs[0] = up->SubVector(0)->SubVector(1)->Data().ptr();
      vecs[1] = up->SubVector(1)->SubVector(1)->Data().ptr();
      surf_db_->WriteVectors(vnames, vecs, true);
        
      vnames[0] = "sub p";
      vnames[1] = "sub T";
      vecs[0] = up->SubVector(0)->SubVector(0)->Data().ptr();
      vecs[1] = up->SubVector(1)->SubVector(0)->Data().ptr();
      domain_db_->WriteVectors(vnames, vecs, true);
    }
  }

  // potentially change for source of water onto ice
  changed |= modify_predictor_for_source_on_ice_(h, up);

  // potentially update faces
  changed |= StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::modify_predictor(h, up);

  return changed;
}


// update the predictor to be physically consistent
bool MPCPermafrost::modify_predictor_for_source_on_ice_(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor, source on ice." << std::endl;

  Teuchos::RCP<const CompositeVector> surf_T = up->SubVector(1)->SubVector(1)->Data();
  const Epetra_MultiVector& surf_T_c = *surf_T->ViewComponent("cell",false);
  const Epetra_MultiVector& mass_src = *S_next_->GetFieldData("surface_mass_source")
      ->ViewComponent("cell",false);
  
  Epetra_MultiVector& surf_p_c = *up->SubVector(0)->SubVector(1)->Data()
      ->ViewComponent("cell",false);
  Epetra_MultiVector& domain_p_f = *up->SubVector(0)->SubVector(0)->Data()
      ->ViewComponent("face",false);
  const Epetra_MultiVector& domain_p_c = *up->SubVector(0)->SubVector(0)->Data()
      ->ViewComponent("cell",false);

  int nchanged_l = 0;
  for (unsigned int sc=0; sc!=surf_T_c.MyLength(); ++sc) {
    if (surf_T_c[0][sc] <= 272.15        // cutoff for 0 water, so surface is all ice
        && mass_src[0][sc] > 0.          // water source is on
        && surf_p_c[0][sc] < 101325.) {  // no ponded water yet
      surf_p_c[0][sc] = 101425.;
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      domain_p_f[0][f] = surf_p_c[0][sc];
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  Modified at surf cell " << sc << ", T = " << surf_T_c[0][sc] << std::endl;
      nchanged_l++;
    } else if (surf_T_c[0][sc] < 273.15
               && surf_p_c[0][sc] > 101325.) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      AmanziMesh::Entity_ID_list cells;
      domain_mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      if (
    }
  }
  int nchanged = nchanged_l;
  domain_mesh_->get_comm()->SumAll(&nchanged_l, &nchanged, 1);
  return nchanged > 0;
}


// updates the preconditioner.  Note this is currently identical to
// MPCCoupledCells::update_precon(), but may change to include surface
// coupling terms.
void MPCPermafrost::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::update_precon(t,up,h);

  if (!decoupled_) {
    S_next_->GetFieldEvaluator(A_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y2_key_);
    S_next_->GetFieldEvaluator(B_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y1_key_);
    Teuchos::RCP<const CompositeVector> dA_dy2 = S_next_->GetFieldData(dA_dy2_key_);
    Teuchos::RCP<const CompositeVector> dB_dy1 = S_next_->GetFieldData(dB_dy1_key_);

    // scale by 1/h
    Teuchos::RCP<Epetra_MultiVector> Ccc =
        Teuchos::rcp(new Epetra_MultiVector(*dA_dy2->ViewComponent("cell",false)));
    (*Ccc) = *dA_dy2->ViewComponent("cell",false);
    Ccc->Scale(1./h);

    Teuchos::RCP<Epetra_MultiVector> Dcc =
        Teuchos::rcp(new Epetra_MultiVector(*dB_dy1->ViewComponent("cell",false)));
    (*Dcc) = *dB_dy1->ViewComponent("cell",false);
    Dcc->Scale(1./h);

    // Assemble the precon, form Schur complement
    mfd_surf_preconditioner_->SetOffDiagonals(Ccc,Dcc);
    mfd_surf_preconditioner_->ComputeSchurComplement();
    mfd_surf_preconditioner_->UpdatePreconditioner();
  }
}


void MPCPermafrost::precon(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

  if (decoupled_) return StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>::precon(u,Pu);

  // make a new TreeVector that is just the subsurface values (by pointer).
  // -- note these const casts are necessary to create the new TreeVector, but
  //    since the TreeVector COULD be const (it is only used in a single method,
  //    in which it is const), const-correctness is not violated here.
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Precon pulling subsurface vectors." << std::endl;

  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector());
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(0)->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(1)->SubVector(0)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector());
  domain_Pu_tv->PushBack(Pu->SubVector(0)->SubVector(0));
  domain_Pu_tv->PushBack(Pu->SubVector(1)->SubVector(0));

  // call the operator's inverse
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Precon applying coupled subsurface operator." << std::endl;
  mfd_surf_preconditioner_->ApplyInverse(*domain_u_tv, domain_Pu_tv.ptr());

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::RCP<const CompositeVector> surf_p = u->SubVector(0)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> surf_Pp = Pu->SubVector(0)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_Pp = Pu->SubVector(0)->SubVector(0)->Data();

    Teuchos::RCP<const CompositeVector> surf_T = u->SubVector(1)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> surf_PT = Pu->SubVector(1)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_PT = Pu->SubVector(1)->SubVector(0)->Data();

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");
    *vo_->os() << " SubSurface precon:" << std::endl;

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(domain_p.ptr());
    vecs.push_back(domain_Pp.ptr());
    vecs.push_back(domain_T.ptr());
    vecs.push_back(domain_PT.ptr());
    domain_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << " Surface precon:" << std::endl;

    vecs[0] = surf_p.ptr();
    vecs[1] = surf_Pp.ptr();
    vecs[2] = surf_T.ptr();
    vecs[3] = surf_PT.ptr();
    surf_db_->WriteVectors(vnames, vecs, true);
  }
  
  // Update source on ice terms
  const Epetra_MultiVector& dWC_dp = *S_next_->GetFieldData("dsurface_water_content_dsurface_pressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& uf = *S_next_->GetFieldData("unfrozen_fraction")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& surf_p_c = *u->SubVector(0)->SubVector(1)->Data()
      ->ViewComponent("cell",false);
  Epetra_MultiVector& domain_Pp_f = *Pu->SubVector(0)->SubVector(0)->Data()
      ->ViewComponent("face",false);
  double dt = S_next_->time() - S_inter_->time();
  for (unsigned int sc=0; sc!=uf.MyLength(); ++sc) {
    if (std::abs(surf_p_c[0][sc]) > 0. && uf[0][sc] == 0.) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      domain_Pp_f[0][f] = surf_p_c[0][sc] / (dWC_dp[0][sc] / dt);
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "Ice precon: dp surface c[" << sc << "] = " << domain_Pp_f[0][f] << " due to res = " << surf_p_c[0][sc] << std::endl;
    }
  }
  
  // Update surface cells and faces on both
  coupled_flow_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(0));
  coupled_flow_pk_->PreconUpdateSurfaceFaces_(u->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(1));
  coupled_energy_pk_->PreconUpdateSurfaceFaces_(u->SubVector(1), Pu->SubVector(1));

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::RCP<const CompositeVector> surf_p = u->SubVector(0)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> surf_Pp = Pu->SubVector(0)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_Pp = Pu->SubVector(0)->SubVector(0)->Data();

    Teuchos::RCP<const CompositeVector> surf_T = u->SubVector(1)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> surf_PT = Pu->SubVector(1)->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_PT = Pu->SubVector(1)->SubVector(0)->Data();

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");
    *vo_->os() << "FINAL Preconditioner application" << std::endl
               << " SubSurface precon:" << std::endl;

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(domain_p.ptr());
    vecs.push_back(domain_Pp.ptr());
    vecs.push_back(domain_T.ptr());
    vecs.push_back(domain_PT.ptr());
    domain_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << " Surface precon:" << std::endl;

    vecs[0] = surf_p.ptr();
    vecs[1] = surf_Pp.ptr();
    vecs[2] = surf_T.ptr();
    vecs[3] = surf_PT.ptr();
    surf_db_->WriteVectors(vnames, vecs, true);
  }

}

bool MPCPermafrost::modify_correction(double h,
        Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Now post-process
  bool flow_modified = coupled_flow_pk_->PreconPostprocess_(res->SubVector(0), Pu->SubVector(0));
  bool energy_modified = coupled_energy_pk_->PreconPostprocess_(res->SubVector(1), Pu->SubVector(1));

  // Update surface cells and faces
  if (flow_modified) {
    coupled_flow_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(0));
    coupled_flow_pk_->PreconUpdateSurfaceFaces_(res->SubVector(0), Pu->SubVector(0));
  }

  if (energy_modified) {
    coupled_energy_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(1));
    coupled_energy_pk_->PreconUpdateSurfaceFaces_(res->SubVector(1), Pu->SubVector(1));
  }
  return flow_modified || energy_modified;
}
} // namespace
