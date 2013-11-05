#include "Teuchos_XMLParameterListHelpers.hpp"

#include "mpc_surface_subsurface_helpers.cc"
#include "permafrost_model.hh"

#include "mpc_permafrost2.hh"

namespace Amanzi {

MPCPermafrost2::MPCPermafrost2(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    StrongMPC<PKPhysicalBDFBase>(plist, FElist, soln) {}


void
MPCPermafrost2::setup(const Teuchos::Ptr<State>& S) {
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");

  // -- turn off PC assembly: <Parameter name="assemble preconditioner" type="bool" value="false"/>
  plist_->sublist("PKs").sublist(names[0]).set("assemble preconditioner", false);
  plist_->sublist("PKs").sublist(names[1]).set("assemble preconditioner", false);
  plist_->sublist("PKs").sublist(names[2]).set("assemble preconditioner", false);
  plist_->sublist("PKs").sublist(names[3]).set("assemble preconditioner", false);

  // -- turn on coupling
  plist_->sublist("PKs").sublist(names[0]).set("coupled to surface via flux", true);
  plist_->sublist("PKs").sublist(names[1]).set("coupled to surface via flux", true);
  plist_->sublist("PKs").sublist(names[2]).set("coupled to subsurface via flux", true);
  plist_->sublist("PKs").sublist(names[3]).set("coupled to subsurface via flux", true);

  // -- set up PC coupling
  plist_->sublist("PKs").sublist(names[0]).sublist("Diffusion PC").set("coupled to surface", true);
  plist_->sublist("PKs").sublist(names[1]).sublist("Diffusion PC").set("coupled to surface", true);
  plist_->sublist("PKs").sublist(names[2]).sublist("Diffusion PC").set("TPFA", true);
  plist_->sublist("PKs").sublist(names[3]).sublist("Diffusion PC").set("TPFA", true);
  Teuchos::writeParameterListToXmlOStream(*plist_, std::cout);

  // grab the meshes
  surf_mesh_ = S->GetMesh("surface");
  domain_mesh_ = S->GetMesh();

  // cast the PKs
  domain_flow_pk_ = sub_pks_[0];
  domain_energy_pk_ = sub_pks_[1];
  surf_flow_pk_ = sub_pks_[2];
  surf_energy_pk_ = sub_pks_[3];

  // call the MPC's setup, which calls the sub-pk's setups
  StrongMPC<PKPhysicalBDFBase>::setup(S);
  
  // require the coupling fields, claim ownership
  S->RequireField("surface_subsurface_flux", name_)
      ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("surface_subsurface_energy_flux", name_)
      ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // create the preconditioner
  Teuchos::ParameterList& pc_sublist = plist_->sublist("Coupled PC");
  precon_ = Teuchos::rcp(new Operators::MatrixMFD_Coupled_Surf(pc_sublist, domain_mesh_));

  pc_flow_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(
      domain_flow_pk_->preconditioner());
  ASSERT(pc_flow_ != Teuchos::null);
  pc_energy_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(
      domain_energy_pk_->preconditioner());
  ASSERT(pc_energy_ != Teuchos::null);
  pc_surf_flow_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(
      surf_flow_pk_->preconditioner());
  ASSERT(pc_surf_flow_ != Teuchos::null);
  pc_surf_energy_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(
      surf_energy_pk_->preconditioner());
  ASSERT(pc_surf_energy_ != Teuchos::null);

  pc_flow_->SetSurfaceOperator(pc_surf_flow_);
  pc_energy_->SetSurfaceOperator(pc_surf_energy_);
  // must re-symbolic assemble surf operators, now that they have a surface operator
  pc_flow_->SymbolicAssembleGlobalMatrices();
  pc_energy_->SymbolicAssembleGlobalMatrices();

  precon_->SetSubBlocks(pc_flow_, pc_energy_);
  precon_->SetSurfaceOperators(pc_surf_flow_, pc_surf_energy_);
  precon_->SymbolicAssembleGlobalMatrices();
  precon_->InitPreconditioner();

  // set up the EWC delegate
  ewc_ = Teuchos::rcp(new MPCDelegateEWC(*plist_));
  Teuchos::RCP<PermafrostModel> model = Teuchos::rcp(new PermafrostModel());
  ewc_->set_model(model);
  ewc_->setup(S);  

  // set up the Water delegate
  water_ = Teuchos::rcp(new MPCDelegateWater(plist_, 0, 2));
  consistent_cells_ =
      plist_->get<bool>("ensure consistent cells after face updates", false);
  
  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  surf_db_ = surf_flow_pk_->debugger();
}

void
MPCPermafrost2::initialize(const Teuchos::Ptr<State>& S) {
  // initialize coupling terms
  S->GetFieldData("surface_subsurface_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_flux", name_)->set_initialized();
  S->GetFieldData("surface_subsurface_energy_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_energy_flux", name_)->set_initialized();


  StrongMPC<PKPhysicalBDFBase>::initialize(S);
  ewc_->initialize(S);
}

void
MPCPermafrost2::set_states(const Teuchos::RCP<const State>& S,
                           const Teuchos::RCP<State>& S_inter,
                           const Teuchos::RCP<State>& S_next) {
  StrongMPC<PKPhysicalBDFBase>::set_states(S,S_inter,S_next);
  ewc_->set_states(S,S_inter,S_next);
}

// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk fun().
void
MPCPermafrost2::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  solution_to_state(u_new, S_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->fun(t_old, t_new, u_old->SubVector(2),
                     u_new->SubVector(2), g->SubVector(2));

  // The residual of the surface flow equation provides the mass flux from
  // subsurface to surface.
  Epetra_MultiVector& source = *S_next_->GetFieldData("surface_subsurface_flux",
          name_)->ViewComponent("cell",false);
  source = *g->SubVector(2)->Data()->ViewComponent("cell",false);
  
  // The exception to this is if the surface unfrozen fraction is 0 and the
  // flux direction is inward, in which case the rel perm will be zero, and no
  // flux is allowed.
  const Epetra_MultiVector& uf_frac = *S_next_->GetFieldData("unfrozen_fraction")
      ->ViewComponent("cell",false);
  for (unsigned int sc=0; sc!=source.MyLength(); ++sc) {
    if (uf_frac[0][sc] == 0. && source[0][sc] < 0.) {
      source[0][sc] = 0.;
    }
  }

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->fun(t_old, t_new, u_old->SubVector(0),
                       u_new->SubVector(0), g->SubVector(0));

  // Clobber the subsurface face's residual if it gets hit with zero rel perm,
  // or else clobber the surface cell's residual as it is taken as a flux into
  // the subsurface.
  Epetra_MultiVector& domain_g_f = *g->SubVector(0)
      ->Data()->ViewComponent("face",false);
  Epetra_MultiVector& surf_g_c = *g->SubVector(2)
      ->Data()->ViewComponent("cell",false);
  for (unsigned int sc=0; sc!=source.MyLength(); ++sc) {
    if (uf_frac[0][sc] == 0. && source[0][sc] < 0.) {
      AmanziMesh::Entity_ID f =
          surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      domain_g_f[0][f] = 0.;
    } else {
      surf_g_c[0][sc] = 0.;
    }
  }
  
  // Now that fluxes are done, do energy.
  // Evaluate the surface energy residual
  surf_energy_pk_->fun(t_old, t_new, u_old->SubVector(3),
                     u_new->SubVector(3), g->SubVector(3));

  // The residual of the surface energy equation provides the diffusive energy
  // flux from subsurface to surface.
  Epetra_MultiVector& esource =
      *S_next_->GetFieldData("surface_subsurface_energy_flux", name_)
      ->ViewComponent("cell",false);
  esource = *g->SubVector(3)->Data()->ViewComponent("cell",false);
  
  // Evaluate the subsurface energy residual.
  domain_energy_pk_->fun(t_old, t_new, u_old->SubVector(1),
                     u_new->SubVector(1), g->SubVector(1));

  // Clobber the surface cell's energy residual, as it is taken as a diffusive
  // flux into the subsurface.
  g->SubVector(3)->PutScalar(0.);

}

// -- Apply preconditioner to u and returns the result in Pu.
void
MPCPermafrost2::precon(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon application:" << std::endl;

  // make a new TreeVector that is just the subsurface values (by pointer).
  // -- note these const casts are necessary to create the new TreeVector, but
  //    since the TreeVector COULD be const (it is only used in a single method,
  //    in which it is const), const-correctness is not violated here.
  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector());
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(u->SubVector(1)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector());
  domain_Pu_tv->PushBack(Pu->SubVector(0));
  domain_Pu_tv->PushBack(Pu->SubVector(1));

  // call the operator's inverse
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying coupled subsurface operator." << std::endl;
  precon_->ApplyInverse(*domain_u_tv, domain_Pu_tv.ptr());

  // Copy subsurface face corrections to surface cell corrections
  CopySubsurfaceToSurface(*Pu->SubVector(0)->Data(),
                          Pu->SubVector(2)->Data().ptr());
  CopySubsurfaceToSurface(*Pu->SubVector(1)->Data(),
                          Pu->SubVector(3)->Data().ptr());

  // Derive surface face corrections.
  UpdateConsistentFaceCorrectionWater_(u, Pu);
  pc_surf_energy_->UpdateConsistentFaceCorrection(*u->SubVector(3)->Data(),
          Pu->SubVector(3)->Data().ptr());

  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> domain_Pp = Pu->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_PT = Pu->SubVector(1)->Data();

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr());
    vecs.push_back(Pu->SubVector(0)->Data().ptr());
    vecs.push_back(u->SubVector(1)->Data().ptr());
    vecs.push_back(Pu->SubVector(1)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = u->SubVector(2)->Data().ptr();
    vecs[1] = Pu->SubVector(2)->Data().ptr();
    vecs[2] = u->SubVector(3)->Data().ptr();
    vecs[3] = Pu->SubVector(3)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);    
  }
}

// -- Update the preconditioner.
void
MPCPermafrost2::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // Create the on-diagonal block PCs
  StrongMPC<PKPhysicalBDFBase>::update_precon(t,up,h);

  // Add the off-diagonal blocks.
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "temperature");
  Teuchos::RCP<const CompositeVector> dWCdT_domain =
      S_next_->GetFieldData("dwater_content_dtemperature");

  S_next_->GetFieldEvaluator("energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "pressure");
  Teuchos::RCP<const CompositeVector> dEdp_domain =
      S_next_->GetFieldData("denergy_dpressure");

  S_next_->GetFieldEvaluator("surface_water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_temperature");
  Teuchos::RCP<const CompositeVector> dWCdT_surf =
      S_next_->GetFieldData("dsurface_water_content_dsurface_temperature");

  S_next_->GetFieldEvaluator("surface_energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_pressure");
  Teuchos::RCP<const CompositeVector> dEdp_surf =
      S_next_->GetFieldData("dsurface_energy_dsurface_pressure");

  precon_->SetOffDiagonals(dWCdT_domain->ViewComponent("cell",false),
                           dEdp_domain->ViewComponent("cell",false),
                           //                           dWCdT_surf->ViewComponent("cell",false),
                           //                           dEdp_surf->ViewComponent("cell",false),
                           Teuchos::null,
                           Teuchos::null,
                           1./h);

  // Assemble the PC
  precon_->ComputeSchurComplement();
  precon_->UpdatePreconditioner();
}

// -- Modify the predictor.
bool
MPCPermafrost2::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  return false;
}

// -- Modify the correction.
bool
MPCPermafrost2::modify_correction(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();
  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "NKA'd, PC'd correction." << std::endl;

    Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> domain_Pp = du->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_PT = du->SubVector(1)->Data();

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());
    vecs.push_back(u->SubVector(1)->Data().ptr());
    vecs.push_back(du->SubVector(1)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = u->SubVector(2)->Data().ptr();
    vecs[1] = du->SubVector(2)->Data().ptr();
    vecs[2] = u->SubVector(3)->Data().ptr();
    vecs[3] = du->SubVector(3)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);    
  }


  bool modified = water_->ModifyCorrection(h, res, u, du);
  if (modified && consistent_cells_) {
    // Derive subsurface cell corrections.
    pc_flow_->UpdateConsistentCellCorrection(
        *u->SubVector(0)->Data(),
        du->SubVector(0)->Data().ptr());
  }

  if (modified) {
    // Copy subsurface face corrections to surface cell corrections
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
                            du->SubVector(2)->Data().ptr());

    // Derive surface face corrections.
    UpdateConsistentFaceCorrectionWater_(res, du);
  }

  // dump to screen
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Modified correction." << std::endl;

    Teuchos::RCP<const CompositeVector> domain_p = u->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> domain_Pp = du->SubVector(0)->Data();
    Teuchos::RCP<const CompositeVector> domain_T = u->SubVector(1)->Data();
    Teuchos::RCP<const CompositeVector> domain_PT = du->SubVector(1)->Data();

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());
    vecs.push_back(u->SubVector(1)->Data().ptr());
    vecs.push_back(du->SubVector(1)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = u->SubVector(2)->Data().ptr();
    vecs[1] = du->SubVector(2)->Data().ptr();
    vecs[2] = u->SubVector(3)->Data().ptr();
    vecs[3] = du->SubVector(3)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);    
  }
  
  return modified;
}

void
MPCPermafrost2::UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
        const Teuchos::RCP<TreeVector>& Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();

  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(2)->Data();
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(2)->Data();

  // Calculate delta h on the surface
  Teuchos::RCP<CompositeVector> surf_Ph = Teuchos::rcp(new CompositeVector(*surf_Pu));
  surf_Ph->PutScalar(0.);

  // old ponded depth
  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
  *surf_Ph->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

  // new ponded depth
  S_next_->GetFieldData("surface_pressure", sub_pks_[2]->name())
      ->ViewComponent("cell",false)->Update(-1., surf_Pu_c, 1.);
  sub_pks_[2]->changed_solution();
  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);

  // put delta ponded depth into surf_Ph_cell
  surf_Ph->ViewComponent("cell",false)
      ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

  // update delta faces
  pc_surf_flow_->UpdateConsistentFaceCorrection(*surf_u, surf_Ph.ptr());
  *surf_Pu->ViewComponent("face",false) = *surf_Ph->ViewComponent("face",false);

  // revert solution so we don't break things
  S_next_->GetFieldData("surface_pressure",sub_pks_[2]->name())
      ->ViewComponent("cell",false)->Update(1., surf_Pu_c, 1.);
  sub_pks_[2]->changed_solution();
}


} // namespace
