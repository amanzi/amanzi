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

#include "mpc_surface_subsurface_flux_coupler.hh"
#include "mpc_permafrost.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

RegisteredPKFactory<MPCPermafrost> MPCPermafrost::reg_("permafrost model");

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void MPCPermafrost::setup(const Teuchos::Ptr<State>& S) {
  // off-diagonal terms needed by MPCCoupledCells
  plist_.set("conserved quantity A", "water_content");
  plist_.set("conserved quantity B", "energy");
  plist_.set("primary variable A", "pressure");
  plist_.set("primary variable B", "temperature");

  plist_.set("mesh key", "domain");
  MPCCoupledCells::setup(S);

  // grab the PKs
  coupled_flow_pk_ =
      Teuchos::rcp_dynamic_cast<MPCSurfaceSubsurfaceFluxCoupler>(sub_pks_[0]);
  ASSERT(coupled_flow_pk_ != Teuchos::null);
  coupled_energy_pk_ =
      Teuchos::rcp_dynamic_cast<MPCSurfaceSubsurfaceFluxCoupler>(sub_pks_[1]);
  ASSERT(coupled_energy_pk_ != Teuchos::null);
}


void MPCPermafrost::precon(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<TreeVector> Pu) {
  if (decoupled_) return StrongMPC::precon(u,Pu);

  // make a new TreeVector that is just the subsurface values (by pointer).
  // -- note these const casts are necessary to create the new TreeVector, but
  //    since the TreeVector COULD be const (it is only used in a single method,
  //    in which it is const), const-correctness is not violated here.
  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector("domain_u_tv"));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(0)->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(1)->SubVector(0)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector("domain_Pu_tv"));
  domain_Pu_tv->PushBack(Pu->SubVector(0)->SubVector(0));
  domain_Pu_tv->PushBack(Pu->SubVector(1)->SubVector(0));

  // call the operator's inverse
  preconditioner_->ApplyInverse(*domain_u_tv, domain_Pu_tv.ptr());

  // Now post-process
  coupled_flow_pk_->PreconPostprocess_(u->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconPostprocess_(u->SubVector(1), Pu->SubVector(1));

  // Update surface cells and faces on both
  coupled_flow_pk_->PreconUpdateSurfaceCells_(u->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceCells_(u->SubVector(1), Pu->SubVector(1));
  coupled_flow_pk_->PreconUpdateSurfaceFaces_(u->SubVector(0), Pu->SubVector(0));
  coupled_energy_pk_->PreconUpdateSurfaceFaces_(u->SubVector(1), Pu->SubVector(1));

#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  Teuchos::RCP<const CompositeVector> surf_p = u->SubVector(0)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> surf_Pp = Pu->SubVector(0)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_Pp = Pu->SubVector(0)->SubVector(0)->data();

  Teuchos::RCP<const CompositeVector> surf_T = u->SubVector(1)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> surf_PT = Pu->SubVector(1)->SubVector(1)->data();
  Teuchos::RCP<const CompositeVector> domain_PT = Pu->SubVector(1)->SubVector(0)->data();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    AmanziMesh::Entity_ID_List fnums1,fnums0;
    std::vector<int> dirs;
    domain_p->mesh()->cell_get_faces_and_dirs(c0_, &fnums0, &dirs);
    domain_p->mesh()->cell_get_faces_and_dirs(c1_, &fnums1, &dirs);

    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;
    *out_ << "  p0: " << (*domain_p)("cell",c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*domain_p)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  p1: " << (*domain_p)("cell",c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*domain_p)("face", fnums1[n]);
    *out_ << std::endl;

    *out_ << "  PC*p0: " << (*domain_Pp)("cell",c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*domain_Pp)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  PC*p1: " << (*domain_Pp)("cell",c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*domain_Pp)("face", fnums1[n]);
    *out_ << std::endl;

    *out_ << "  ---" << std::endl;

    *out_ << "  T0: " << (*domain_T)("cell",c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*domain_T)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  T1: " << (*domain_T)("cell",c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*domain_T)("face", fnums1[n]);
    *out_ << std::endl;

    *out_ << "  PC*T0: " << (*domain_PT)("cell",c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*domain_PT)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  PC*T1: " << (*domain_PT)("cell",c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*domain_PT)("face", fnums1[n]);
    *out_ << std::endl;
    *out_ << "  ---" << std::endl;

  }

  if (coupled_flow_pk_->surf_c0_ < surf_p->size("cell",false) && coupled_flow_pk_->surf_c1_ < surf_p->size("cell",false)) {
    //  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
     AmanziMesh::Entity_ID_List fnums1,fnums0;
     std::vector<int> dirs;
     surf_p->mesh()->cell_get_faces_and_dirs(coupled_flow_pk_->surf_c0_, &fnums0, &dirs);
     surf_p->mesh()->cell_get_faces_and_dirs(coupled_flow_pk_->surf_c1_, &fnums1, &dirs);

    *out_ << " Surface precon:" << std::endl;
    *out_ << "  p0: " << (*surf_p)("cell",coupled_flow_pk_->surf_c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*surf_p)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  p1: " << (*surf_p)("cell",coupled_flow_pk_->surf_c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*surf_p)("face", fnums1[n]);
    *out_ << std::endl;

    *out_ << "  PC*p0: " << (*surf_Pp)("cell",coupled_flow_pk_->surf_c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*surf_Pp)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  PC*p1: " << (*surf_Pp)("cell",coupled_flow_pk_->surf_c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*surf_Pp)("face", fnums1[n]);
    *out_ << std::endl;

    *out_ << "  ---" << std::endl;

    *out_ << "  T0: " << (*surf_T)("cell",coupled_flow_pk_->surf_c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*surf_T)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  T1: " << (*surf_T)("cell",coupled_flow_pk_->surf_c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*surf_T)("face", fnums1[n]);
    *out_ << std::endl;

    *out_ << "  PC*T0: " << (*surf_PT)("cell",coupled_flow_pk_->surf_c0_);
    for (int n=0;n!=fnums0.size();++n)
      *out_ << ", " << (*surf_PT)("face", fnums0[n]);
    *out_ << std::endl;

    *out_ << "  PC*T1: " << (*surf_PT)("cell",coupled_flow_pk_->surf_c1_);
    for (int n=0;n!=fnums1.size();++n)
      *out_ << ", " << (*surf_PT)("face", fnums1[n]);
    *out_ << std::endl;
    *out_ << "  ---" << std::endl;
  }
#endif

}

} // namespace
