#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "mpc_delegate_ewc_surface.hh"
#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "permafrost_model.hh"
#include "surface_ice_model.hh"

#include "mpc_permafrost3.hh"

namespace Amanzi {

MPCPermafrost3::MPCPermafrost3(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    StrongMPC<PKPhysicalBDFBase>(plist, FElist, soln) {}


void
MPCPermafrost3::setup(const Teuchos::Ptr<State>& S) {
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

  // Create the preconditioner.
  // -- Allocate matrix.
  Teuchos::ParameterList& pc_sublist = plist_->sublist("Coupled PC");
  precon_ = Teuchos::rcp(new Operators::MatrixMFD_Coupled_Surf(pc_sublist, domain_mesh_));

  // -- Collect the sub-blocks.
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

  // -- Subsurface blocks include their surface operators.
  pc_flow_->SetSurfaceOperator(pc_surf_flow_);
  pc_energy_->SetSurfaceOperator(pc_surf_energy_);

  // -- must re-symbolic assemble surf operators, now that they have a surface operator
  pc_flow_->SymbolicAssembleGlobalMatrices();
  pc_energy_->SymbolicAssembleGlobalMatrices();

  // -- finally, set the sub blocks in the coupled PC
  precon_->SetSubBlocks(pc_flow_, pc_energy_);
  precon_->SetSurfaceOperators(pc_surf_flow_, pc_surf_energy_);
  precon_->SymbolicAssembleGlobalMatrices();
  precon_->InitPreconditioner();

  // Potential create the solver
  // NOTE: for this to be enabled, we must first implement the forward
  // operator's assembly for MatrixMFD_Coupled_Surf -- etc
  // if (plist_->isSublist("Coupled Solver")) {
  //   Teuchos::ParameterList linsolve_sublist = plist_->sublist("Coupled Solver");
  //   AmanziSolvers::LinearOperatorFactory<TreeMatrix,TreeVector,TreeVectorSpace> fac;
  //   lin_solver_ = fac.Create(linsolve_sublist, precon_);
  // } else {
    lin_solver_ = precon_;
  // }
  
  // set up the EWC delegates
  Teuchos::RCP<Teuchos::ParameterList> sub_ewc_list = Teuchos::sublist(plist_, "subsurface ewc delegate");
  sub_ewc_list->set("PK name", name_);
  sub_ewc_list->set("domain key", "");
  sub_ewc_ = Teuchos::rcp(new MPCDelegateEWCSubsurface(*sub_ewc_list));
  Teuchos::RCP<PermafrostModel> sub_model = Teuchos::rcp(new PermafrostModel());
  sub_ewc_->set_model(sub_model);
  sub_ewc_->setup(S);  

  Teuchos::RCP<Teuchos::ParameterList> surf_ewc_list = Teuchos::sublist(plist_, "surface ewc delegate");
  surf_ewc_list->set("PK name", name_);
  surf_ewc_list->set("domain key", "surface");
  surf_ewc_ = Teuchos::rcp(new MPCDelegateEWCSurface(*surf_ewc_list));
  Teuchos::RCP<SurfaceIceModel> surf_model = Teuchos::rcp(new SurfaceIceModel());
  surf_ewc_->set_model(surf_model);
  surf_ewc_->setup(S);  
  
  // set up the Water delegate
  Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
  water_ = Teuchos::rcp(new MPCDelegateWater(water_list));
  water_->set_indices(0,2,1,3);

  // set up our own predictors
  consistent_cells_ =
      plist_->get<bool>("ensure consistent cells after face updates", false);
  
  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  surf_db_ = surf_flow_pk_->debugger();
}

void
MPCPermafrost3::initialize(const Teuchos::Ptr<State>& S) {
  // initialize coupling terms
  S->GetFieldData("surface_subsurface_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_flux", name_)->set_initialized();
  S->GetFieldData("surface_subsurface_energy_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_energy_flux", name_)->set_initialized();

  // Initialize all sub PKs.
  MPC<PKPhysicalBDFBase>::initialize(S);

  // ensure continuity of ICs... surface takes precedence.
  CopySurfaceToSubsurface(*S->GetFieldData("surface_pressure", sub_pks_[2]->name()),
                          S->GetFieldData("pressure", sub_pks_[0]->name()).ptr());
  CopySurfaceToSubsurface(*S->GetFieldData("surface_temperature", sub_pks_[3]->name()),
                          S->GetFieldData("temperature", sub_pks_[1]->name()).ptr());

  // initialize delegates
  surf_ewc_->initialize(S);
  sub_ewc_->initialize(S);

  // Initialize my timestepper.
  PKBDFBase::initialize(S);
}

void
MPCPermafrost3::set_states(const Teuchos::RCP<const State>& S,
                           const Teuchos::RCP<State>& S_inter,
                           const Teuchos::RCP<State>& S_next) {
  StrongMPC<PKPhysicalBDFBase>::set_states(S,S_inter,S_next);
  surf_ewc_->set_states(S,S_inter,S_next);
  sub_ewc_->set_states(S,S_inter,S_next);
  water_->set_states(S,S_inter,S_next);
}


void
MPCPermafrost3::commit_state(double dt, const Teuchos::RCP<State>& S) {
  StrongMPC<PKPhysicalBDFBase>::commit_state(dt,S);
  if (surf_ewc_ != Teuchos::null) surf_ewc_->commit_state(dt,S);
  if (sub_ewc_ != Teuchos::null) sub_ewc_->commit_state(dt,S);
}


// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk fun().
void
MPCPermafrost3::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
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

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->fun(t_old, t_new, u_old->SubVector(0),
                       u_new->SubVector(0), g->SubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(2)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  // Now that mass fluxes are done, do energy.
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

  // All energy fluxes have been taken by the subsurface.
  g->SubVector(3)->Data()->ViewComponent("cell",false)->PutScalar(0.);
}

// -- Apply preconditioner to u and returns the result in Pu.
void
MPCPermafrost3::precon(Teuchos::RCP<const TreeVector> u,
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
  lin_solver_->ApplyInverse(*domain_u_tv, *domain_Pu_tv);

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
MPCPermafrost3::update_precon(double t,
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

  // ALWAYS 0!
  // S_next_->GetFieldEvaluator("surface_water_content")
  //     ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_temperature");
  // Teuchos::RCP<const CompositeVector> dWCdT_surf =
  //     S_next_->GetFieldData("dsurface_water_content_dsurface_temperature");

  S_next_->GetFieldEvaluator("surface_energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_pressure");
  Teuchos::RCP<const CompositeVector> dEdp_surf =
      S_next_->GetFieldData("dsurface_energy_dsurface_pressure");

  precon_->SetOffDiagonals(dWCdT_domain->ViewComponent("cell",false),
                           dEdp_domain->ViewComponent("cell",false),
                           Teuchos::null, // dWC_dT = 0
                           dEdp_surf->ViewComponent("cell",false),
                           1./h);

  // Assemble the PC
  precon_->ComputeSchurComplement();

  // // dump the schur complement
  // Teuchos::RCP<const Epetra_FEVbrMatrix> sc = precon_->Schur();
  // std::stringstream filename_s;
  // filename_s << "schur_" << S_next_->cycle() << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);

  precon_->UpdatePreconditioner();
}

// -- Modify the predictor.
bool
MPCPermafrost3::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  bool modified = false;

  // Surface EWC
  // -- make a new TreeVector that is just the surface values (by pointer).
  Teuchos::RCP<TreeVector> surf_u = Teuchos::rcp(new TreeVector());
  surf_u->PushBack(u->SubVector(2));
  surf_u->PushBack(u->SubVector(3));
  modified |= surf_ewc_->modify_predictor(h,surf_u);
  // -- copy surf --> sub
  if (modified) {
    CopySurfaceToSubsurface(*u->SubVector(2)->Data(), u->SubVector(0)->Data().ptr());
    CopySurfaceToSubsurface(*u->SubVector(3)->Data(), u->SubVector(1)->Data().ptr());
  }

  // Subsurface EWC
  // -- make a new TreeVector that is just the subsurface values (by pointer).
  Teuchos::RCP<TreeVector> sub_u = Teuchos::rcp(new TreeVector());
  sub_u->PushBack(u->SubVector(0));
  sub_u->PushBack(u->SubVector(1));
  modified |= sub_ewc_->modify_predictor(h,sub_u);

  modified |= water_->ModifyPredictor_Heuristic(h, u);
  modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
  modified |= water_->ModifyPredictor_TempFromSource(h, u);

  modified |= StrongMPC<PKPhysicalBDFBase>::modify_predictor(h, u);
  return modified;
}

// -- Modify the correction.
bool
MPCPermafrost3::modify_correction(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();
  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "NKA'd, PC'd correction." << std::endl;

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(res->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());
    vecs.push_back(res->SubVector(1)->Data().ptr());
    vecs.push_back(du->SubVector(1)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->SubVector(2)->Data().ptr();
    vecs[1] = du->SubVector(2)->Data().ptr();
    vecs[2] = res->SubVector(3)->Data().ptr();
    vecs[3] = du->SubVector(3)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);    
  }


  // modify correction using water approaches
  int n_modified = 0;
  n_modified += water_->ModifyCorrection_WaterFaceLimiter(h, res, u, du);
  double damping = water_->ModifyCorrection_WaterSpurtDamp(h, res, u, du);
  n_modified += water_->ModifyCorrection_WaterSpurtCap(h, res, u, du, damping);

  // -- accumulate globally
  int n_modified_l = n_modified;
  u->SubVector(0)->Data()->Comm().SumAll(&n_modified_l, &n_modified, 1);
  bool modified = (n_modified > 0) || (damping < 1.);

  // -- calculate consistent subsurface cells
  if (modified) {
    if (consistent_cells_) {
      // Derive subsurface cell corrections.
      pc_flow_->UpdateConsistentCellCorrection(
          *u->SubVector(0)->Data(),
          du->SubVector(0)->Data().ptr());
    }

    // Copy subsurface face corrections to surface cell corrections
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
                            du->SubVector(2)->Data().ptr());
  }

  // -- damping can be then applied to energy corrections well
  // if (damping < 1.) {
  //   du->SubVector(1)->Scale(damping);
  //   du->SubVector(3)->Scale(damping);
  // }

  // modify correction for dumping water onto a frozen surface
  //  n_modified = ModifyCorrection_FrozenSurface_(h, res, u, du);
  // -- accumulate globally
  // n_modified_l = n_modified;
  // u->SubVector(0)->Data()->Comm().SumAll(&n_modified_l, &n_modified, 1);
  // modified |= (n_modified > 0) || (damping < 1.);

  if (modified) {
    // Derive surface face corrections.
    UpdateConsistentFaceCorrectionWater_(res, du);
  }

  // dump to screen
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Modified correction." << std::endl;

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");
    vnames.push_back("T");
    vnames.push_back("PC*T");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(res->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());
    vecs.push_back(res->SubVector(1)->Data().ptr());
    vecs.push_back(du->SubVector(1)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->SubVector(2)->Data().ptr();
    vecs[1] = du->SubVector(2)->Data().ptr();
    vecs[2] = res->SubVector(3)->Data().ptr();
    vecs[3] = du->SubVector(3)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);    
  }
  
  return modified;
}

void
MPCPermafrost3::UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
        const Teuchos::RCP<TreeVector>& Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();

  Teuchos::RCP<CompositeVector> surf_Pp = Pu->SubVector(2)->Data();
  Epetra_MultiVector& surf_Pp_c = *surf_Pp->ViewComponent("cell",false);

  Teuchos::RCP<CompositeVector> surf_PT = Pu->SubVector(3)->Data();
  Epetra_MultiVector& surf_PT_c = *surf_PT->ViewComponent("cell",false);

  Teuchos::RCP<const CompositeVector> surf_p = u->SubVector(2)->Data();

  // Calculate delta h on the surface
  Teuchos::RCP<CompositeVector> surf_Ph = Teuchos::rcp(new CompositeVector(*surf_Pp));
  surf_Ph->PutScalar(0.);

  // old ponded depth
  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
  *surf_Ph->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

  // new ponded depth
  Teuchos::RCP<TreeVector> tv_p = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<CompositeVector> cv_p = S_next_->GetFieldData("surface_pressure", sub_pks_[2]->name());
  cv_p->ViewComponent("cell",false)->Update(-1., surf_Pp_c, 1.);
  tv_p->SetData(cv_p);

  Teuchos::RCP<TreeVector> tv_T = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<CompositeVector> cv_T = S_next_->GetFieldData("surface_temperature", sub_pks_[3]->name());
  cv_T->ViewComponent("cell",false)->Update(-1., surf_PT_c, 1.);
  tv_T->SetData(cv_T);

  sub_pks_[2]->changed_solution();
  sub_pks_[3]->changed_solution();

  if (sub_pks_[2]->is_admissible(tv_p) &&
      sub_pks_[3]->is_admissible(tv_T)) {
    S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);

    // put delta ponded depth into surf_Ph_cell
    surf_Ph->ViewComponent("cell",false)
        ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

    // update delta faces
    pc_surf_flow_->UpdateConsistentFaceCorrection(*surf_p, surf_Ph.ptr());
    *surf_Pp->ViewComponent("face",false) = *surf_Ph->ViewComponent("face",false);
  }

  // revert solution so we don't break things
  S_next_->GetFieldData("surface_pressure",sub_pks_[2]->name())
      ->ViewComponent("cell",false)->Update(1., surf_Pp_c, 1.);
  sub_pks_[2]->changed_solution();
  S_next_->GetFieldData("surface_temperature",sub_pks_[3]->name())
      ->ViewComponent("cell",false)->Update(1., surf_PT_c, 1.);
  sub_pks_[3]->changed_solution();
}


int
MPCPermafrost3::ModifyCorrection_FrozenSurface_(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  const Epetra_MultiVector& res_surf_c = *res->SubVector(2)->Data()->ViewComponent("cell",false);
  const Epetra_MultiVector& u_surf_c = *u->SubVector(2)->Data()->ViewComponent("cell",false);

  Epetra_MultiVector& du_surf_c = *du->SubVector(2)->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& du_domain_f = *du->SubVector(0)->Data()->ViewComponent("face",false);

  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  S_next_->GetFieldEvaluator("surface_water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_pressure");
  const Epetra_MultiVector& dWC_dp_surf_c =
      *S_next_->GetFieldData("dsurface_water_content_dsurface_pressure")
      ->ViewComponent("cell",false);

  int nmodified = 0;
  for (unsigned int sc=0; sc!=du_surf_c.MyLength(); ++sc) {
    if (std::abs(res_surf_c[0][sc]) > 0.) {
      nmodified++;

      du_surf_c[0][sc] = u_surf_c[0][sc] - (patm - res_surf_c[0][sc] / (dWC_dp_surf_c[0][sc] / h));
      AmanziMesh::Entity_ID f =
          surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      du_domain_f[0][f] = du_surf_c[0][sc];
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  Source on frozen surface: res = " << res_surf_c[0][sc] << ", dWC_dp = " << dWC_dp_surf_c[0][sc]
                   << ", du = " << du_surf_c[0][sc] << std::endl;
    }
  }
  return nmodified;
}


} // namespace
