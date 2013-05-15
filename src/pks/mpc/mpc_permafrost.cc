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

RegisteredPKFactory<MPCPermafrost> MPCPermafrost::reg_("permafrost model");

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void MPCPermafrost::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC::setup(S);

  // option to decoupled and revert to StrongMPC
  decoupled_ = plist_.get<bool>("decoupled",false);

  // off-diagonal terms needed by MPCCoupledCells
  A_key_ = plist_.get<std::string>("conserved quantity A", "water_content");
  B_key_ = plist_.get<std::string>("conserved quantity B", "energy");
  y1_key_ = plist_.get<std::string>("primary variable A", "pressure");
  y2_key_ = plist_.get<std::string>("primary variable B", "temperature");
  dA_dy2_key_ = std::string("d")+A_key_+std::string("_d")+y2_key_;
  dB_dy1_key_ = std::string("d")+B_key_+std::string("_d")+y1_key_;

  Key mesh_key = plist_.get<std::string>("mesh key", "domain");
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(mesh_key);

  Key surf_mesh_key = plist_.get<std::string>("surface mesh key", "surface");
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh(surf_mesh_key);

  Teuchos::ParameterList pc_sublist = plist_.sublist("Coupled PC");
  mfd_surf_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Coupled_Surf(
      pc_sublist, mesh, surf_mesh));
  preconditioner_ = mfd_surf_preconditioner_;

  // Set the subblocks.  Note these are the flux-coupled PCs, which are
  // MatrixMFD_Surfs.
  Teuchos::RCP<Operators::Matrix> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Matrix> pcB = sub_pks_[1]->preconditioner();

#ifdef ENABLE_DBC
  Teuchos::RCP<Operators::MatrixMFD_Surf> pcB_mfd =
      Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(pcB);
  ASSERT(pcB_mfd != Teuchos::null);
  Teuchos::RCP<Operators::MatrixMFD_Surf> pcA_mfd =
      Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(pcA);
  ASSERT(pcA_mfd != Teuchos::null);
#else
  Teuchos::RCP<Operators::MatrixMFD_Surf> pcA_mfd =
      Teuchos::rcp_static_cast<Operators::MatrixMFD_Surf>(pcA);
  Teuchos::RCP<Operators::MatrixMFD_Surf> pcB_mfd =
      Teuchos::rcp_static_cast<Operators::MatrixMFD_Surf>(pcB);
#endif

  Teuchos::RCP<Operators::MatrixMFD_TPFA> pcA_surf;
  pcA_mfd->GetSurfaceOperator(pcA_surf);
  Teuchos::RCP<Operators::MatrixMFD_TPFA> pcB_surf;
  pcB_mfd->GetSurfaceOperator(pcB_surf);

  mfd_surf_preconditioner_->SetSubBlocks(pcA_mfd, pcB_mfd);
  mfd_surf_preconditioner_->SetSurfaceOperators(pcA_surf, pcB_surf);

  // setup and initialize the preconditioner
  mfd_surf_preconditioner_->SymbolicAssembleGlobalMatrices();
  mfd_surf_preconditioner_->InitPreconditioner();

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_.get<std::string>("predictor type", "none");
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
    ewc_ = Teuchos::rcp(new MPCDelegateEWC(plist_));

    Teuchos::RCP<PermafrostModel> model = Teuchos::rcp(new PermafrostModel());
    ewc_->set_model(model);
    ewc_->setup(S);
  }

  // grab the PKs
  coupled_flow_pk_ =
      Teuchos::rcp_dynamic_cast<MPCSurfaceSubsurfaceFluxCoupler>(sub_pks_[0]);
  ASSERT(coupled_flow_pk_ != Teuchos::null);
  coupled_energy_pk_ =
      Teuchos::rcp_dynamic_cast<MPCSurfaceSubsurfaceFluxCoupler>(sub_pks_[1]);
  ASSERT(coupled_energy_pk_ != Teuchos::null);
}


void MPCPermafrost::initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC::initialize(S);
  if (ewc_ != Teuchos::null) ewc_->initialize(S);
}


void MPCPermafrost::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  StrongMPC::set_states(S,S_inter,S_next);
  if (ewc_ != Teuchos::null) ewc_->set_states(S,S_inter,S_next);
}


void MPCPermafrost::commit_state(double dt, const Teuchos::RCP<State>& S) {
  StrongMPC::commit_state(dt,S);
  if (ewc_ != Teuchos::null) ewc_->commit_state(dt,S);
}


// update the predictor to be physically consistent
bool MPCPermafrost::modify_predictor(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = getOSTab();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "Modifying predictor." << std::endl;

    for (std::vector<int>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      AmanziMesh::Entity_ID_List fnums0;
      std::vector<int> dirs;
      up->SubVector(0)->SubVector(0)->data()->mesh()->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

      *out_ << "  old vals: p = " << (*up->SubVector(0)->SubVector(0)->data())("cell",*c0)
            << ", " << (*up->SubVector(0)->SubVector(0)->data())("face",fnums0[0]) << std::endl;
      *out_ << "            T = " << (*up->SubVector(1)->SubVector(0)->data())("cell",*c0)
            << ", " << (*up->SubVector(1)->SubVector(0)->data())("face",fnums0[0]) << std::endl;
    }
  }


  bool changed(false);
  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    // make a new TreeVector that is just the subsurface values (by pointer).
    Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector("domain_u_tv"));
    domain_u_tv->PushBack(up->SubVector(0)->SubVector(0));
    domain_u_tv->PushBack(up->SubVector(1)->SubVector(0));
    changed = ewc_->modify_predictor(h, domain_u_tv);
  }

  // potentially update faces
  changed |= StrongMPC::modify_predictor(h, up);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    for (std::vector<int>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      AmanziMesh::Entity_ID_List fnums0;
      std::vector<int> dirs;
      up->SubVector(0)->SubVector(0)->data()->mesh()->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

      *out_ << "  new vals: p = " << (*up->SubVector(0)->SubVector(0)->data())("cell",*c0)
            << ", " << (*up->SubVector(0)->SubVector(0)->data())("face",fnums0[0]) << std::endl;
      *out_ << "            T = " << (*up->SubVector(1)->SubVector(0)->data())("cell",*c0)
            << ", " << (*up->SubVector(1)->SubVector(0)->data())("face",fnums0[0]) << std::endl;
    }
  }

  return changed;
}


// updates the preconditioner.  Note this is currently identical to
// MPCCoupledCells::update_precon(), but may change to include surface
// coupling terms.
void MPCPermafrost::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
    *out_ << "Precon update at t = " << t << std::endl;

  StrongMPC::update_precon(t,up,h);

  if (!decoupled_) {
    S_next_->GetFieldEvaluator(A_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y2_key_);
    S_next_->GetFieldEvaluator(B_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y1_key_);
    Teuchos::RCP<const CompositeVector> dA_dy2 = S_next_->GetFieldData(dA_dy2_key_);
    Teuchos::RCP<const CompositeVector> dB_dy1 = S_next_->GetFieldData(dB_dy1_key_);

    // scale by 1/h
    Epetra_MultiVector Ccc(*dA_dy2->ViewComponent("cell",false));
    Ccc = *dA_dy2->ViewComponent("cell",false);
    Ccc.Scale(1./h);

    Epetra_MultiVector Dcc(*dB_dy1->ViewComponent("cell",false));
    Dcc = *dB_dy1->ViewComponent("cell",false);
    Dcc.Scale(1./h);

    // Assemble the precon, form Schur complement
    mfd_surf_preconditioner_->ComputeSchurComplement(Ccc, Dcc);
    mfd_surf_preconditioner_->UpdatePreconditioner();
  }
}


void MPCPermafrost::precon(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
    *out_ << "Precon application:" << std::endl;

  if (decoupled_) return StrongMPC::precon(u,Pu);

  // make a new TreeVector that is just the subsurface values (by pointer).
  // -- note these const casts are necessary to create the new TreeVector, but
  //    since the TreeVector COULD be const (it is only used in a single method,
  //    in which it is const), const-correctness is not violated here.
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
    *out_ << "  Precon pulling subsurface vectors." << std::endl;

  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector("domain_u_tv"));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(0)->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(1)->SubVector(0)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector("domain_Pu_tv"));
  domain_Pu_tv->PushBack(Pu->SubVector(0)->SubVector(0));
  domain_Pu_tv->PushBack(Pu->SubVector(1)->SubVector(0));

  // call the operator's inverse
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
    *out_ << "  Precon applying coupled subsurface operator." << std::endl;
  mfd_surf_preconditioner_->ApplyInverse(*domain_u_tv, domain_Pu_tv.ptr());

  // Update surface cells and faces on both
  coupled_flow_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(1));
  coupled_flow_pk_->PreconUpdateSurfaceFaces_(u->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceFaces_(u->SubVector(1), Pu->SubVector(1));

#if DEBUG_FLAG
  Teuchos::RCP<const CompositeVector> surf_p = u->SubVector(0)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> surf_Pp = Pu->SubVector(0)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_Pp = Pu->SubVector(0)->SubVector(0)->data();

  Teuchos::RCP<const CompositeVector> surf_T = u->SubVector(1)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> surf_PT = Pu->SubVector(1)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_PT = Pu->SubVector(1)->SubVector(0)->data();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;

    for (std::vector<int>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      AmanziMesh::Entity_ID_List fnums0;
      std::vector<int> dirs;
      domain_p->mesh()->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

      *out_ << "  p(" << *c0 << "): " << (*domain_p)("cell",*c0);
      for (int n=0;n!=fnums0.size();++n)
        *out_ << ", " << (*domain_p)("face", fnums0[n]);
      *out_ << std::endl;

      *out_ << "  PC*p(" << *c0 << "): " << (*domain_Pp)("cell",*c0);
      for (int n=0;n!=fnums0.size();++n)
        *out_ << ", " << (*domain_Pp)("face", fnums0[n]);
      *out_ << std::endl;

      *out_ << "  ---" << std::endl;

      *out_ << "  T(" << *c0 << "): " << (*domain_T)("cell",*c0);
      for (int n=0;n!=fnums0.size();++n)
        *out_ << ", " << (*domain_T)("face", fnums0[n]);
      *out_ << std::endl;

      *out_ << "  PC*T(" << *c0 << "): " << (*domain_PT)("cell",*c0);
      for (int n=0;n!=fnums0.size();++n)
        *out_ << ", " << (*domain_PT)("face", fnums0[n]);
      *out_ << std::endl;

      *out_ << "  ---" << std::endl;
    }
  }

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << " Surface precon:" << std::endl;

    for (std::vector<int>::const_iterator c0=coupled_flow_pk_->surf_dc_.begin();
         c0!=coupled_flow_pk_->surf_dc_.end(); ++c0) {
      if (*c0 < surf_p->size("cell",false)) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        surf_p->mesh()->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  p(" << *c0 << "): " << (*surf_p)("cell",*c0);
        for (int n=0;n!=fnums0.size();++n)
          *out_ << ", " << (*surf_p)("face", fnums0[n]);
        *out_ << std::endl;

        *out_ << "  PC*p(" << *c0 << "): " << (*surf_Pp)("cell",*c0);
        for (int n=0;n!=fnums0.size();++n)
          *out_ << ", " << (*surf_Pp)("face", fnums0[n]);
        *out_ << std::endl;

        *out_ << "  ---" << std::endl;

        *out_ << "  T(" << *c0 << "): " << (*surf_T)("cell",*c0);
        for (int n=0;n!=fnums0.size();++n)
          *out_ << ", " << (*surf_T)("face", fnums0[n]);
        *out_ << std::endl;

        *out_ << "  PC*T(" << *c0 << "): " << (*surf_PT)("cell",*c0);
        for (int n=0;n!=fnums0.size();++n)
          *out_ << ", " << (*surf_PT)("face", fnums0[n]);
        *out_ << std::endl;

        *out_ << "  ---" << std::endl;
      }
    }
  }
#endif

}

bool MPCPermafrost::modify_correction(double h,
        Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Now post-process
  coupled_flow_pk_->PreconPostprocess_(res->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconPostprocess_(res->SubVector(1), Pu->SubVector(1));

  // Update surface cells and faces on both
  coupled_flow_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceCells_(Pu->SubVector(1));
  coupled_flow_pk_->PreconUpdateSurfaceFaces_(res->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceFaces_(res->SubVector(1), Pu->SubVector(1));
  return true;
}
} // namespace
