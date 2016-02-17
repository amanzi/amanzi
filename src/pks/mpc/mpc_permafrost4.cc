#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"
#include "mpc_delegate_ewc_surface.hh"
#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "permafrost_model.hh"
#include "surface_ice_model.hh"
#include "energy_base.hh"
#include "advection.hh"

#include "mpc_permafrost4.hh"

namespace Amanzi {

MPCPermafrost4::MPCPermafrost4(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    MPCSubsurface(plist, FElist, soln) {}


void
MPCPermafrost4::setup(const Teuchos::Ptr<State>& S) {
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");

  // -- turn on coupling
  plist_->sublist("PKs").sublist(names[0]).set("coupled to surface via flux", true);
  plist_->sublist("PKs").sublist(names[1]).set("coupled to surface via flux", true);
  plist_->sublist("PKs").sublist(names[2]).set("coupled to subsurface via flux", true);
  plist_->sublist("PKs").sublist(names[3]).set("coupled to subsurface via flux", true);

  // -- ensure local ops are suface ops
  plist_->sublist("PKs").sublist(names[2]).sublist("Diffusion PC").set("surface operator", true);
  plist_->sublist("PKs").sublist(names[2]).sublist("Accumulation PC").set("surface operator", true);
  plist_->sublist("PKs").sublist(names[3]).sublist("Diffusion PC").set("surface operator", true);
  plist_->sublist("PKs").sublist(names[3]).sublist("Advection PC").set("surface operator", true);
  plist_->sublist("PKs").sublist(names[3]).sublist("Accumulation PC").set("surface operator", true);
  
  // grab the meshes
  surf_mesh_ = S->GetMesh("surface");
  domain_mesh_ = S->GetMesh();

  // alias the PKs for easier reference
  domain_flow_pk_ = sub_pks_[0];
  domain_energy_pk_ = sub_pks_[1];
  surf_flow_pk_ = sub_pks_[2];
  surf_energy_pk_ = sub_pks_[3];

  // Create the dE_dp block, which will at least have a CELL-based diagonal
  // entry (from subsurface dE/dp) and a FACE-based diagonal entry (from
  // surface dE/dp), but the subsurface will likely create a CELL-only matrix.
  // This can get removed/fixed once there is a better way of
  // creating/amalgamating ops into a single global operator.
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(domain_mesh_)->SetGhosted()
      ->AddComponent("face", AmanziMesh::FACE, 1)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList plist;
  dE_dp_block_ = Teuchos::rcp(new Operators::Operator_FaceCell(cvs, plist));
  
  // call the subsurface setup, which calls the sub-pk's setups and sets up
  // the subsurface block operator
  MPCSubsurface::setup(S);

  // require the coupling fields, claim ownership
  S->RequireField("surface_subsurface_flux", name_)
      ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("surface_subsurface_energy_flux", name_)
      ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  if (precon_type_ != PRECON_NONE) {  
    // Add the (diagonal) surface blocks into the subsurface blocks.

    // For now we have just the basics, but this could get as complex as
    // MPCSubsurface with offdiagonal terms for surface advection, derivatives
    // of surface conductivity with respect to temperature, etc.

    // -- surface flow
    Teuchos::RCP<Operators::Operator> surf_flow_pc = surf_flow_pk_->preconditioner();
    Teuchos::RCP<Operators::Operator> domain_flow_pc = domain_flow_pk_->preconditioner();
    for (Operators::Operator::op_iterator op = surf_flow_pc->OpBegin();
         op != surf_flow_pc->OpEnd(); ++op) {
      domain_flow_pc->OpPushBack(*op);
    }

    // -- surface energy
    Teuchos::RCP<Operators::Operator> surf_energy_pc = surf_energy_pk_->preconditioner();
    Teuchos::RCP<Operators::Operator> domain_energy_pc = domain_energy_pk_->preconditioner();
    for (Operators::Operator::op_iterator op = surf_energy_pc->OpBegin();
         op != surf_energy_pc->OpEnd(); ++op) {
      domain_energy_pc->OpPushBack(*op);
    }

    if (precon_type_ != PRECON_BLOCK_DIAGONAL) {

      // Add off-diagonal blocks for the surface
      // -- derivatives of surface water content with respect to surface temperature
      // -- ALWAYS ZERO!
      // ASSERT(dWC_dT_block_ != Teuchos::null);
      // Teuchos::ParameterList dWC_dT_plist;
      // dWC_dT_plist.set("surface operator", true);
      // dWC_dT_plist.set("entity kind", "cell");
      // dWC_dT_surf_ = Teuchos::rcp(new Operators::OperatorAccumulation(dWC_dT_plist, dWC_dT_block_));

      // -- derivatives of surface energy with respect to surface pressure
      //    For the Operator, we have to create one with the surface mesh,
      //    then push the op into the full (subsurface) operator.
      ASSERT(dE_dp_block_ != Teuchos::null);
      Teuchos::ParameterList dE_dp_plist;
      dE_dp_plist.set("surface operator", true);
      dE_dp_plist.set("entity kind", "cell");
      dE_dp_surf_ = Teuchos::rcp(new Operators::OperatorAccumulation(dE_dp_plist, surf_mesh_));
      dE_dp_block_->OpPushBack(dE_dp_surf_->local_matrices());
    }

    // must now re-symbolic assemble the matrix to get the updated surface parts
    preconditioner_->SymbolicAssembleMatrix();

  }
      
  // set up the Water delegate
  Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
  water_ = Teuchos::rcp(new MPCDelegateWater(water_list));
  water_->set_indices(0,2,1,3);

  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  surf_db_ = surf_flow_pk_->debugger();

  water_->set_db(surf_db_);
}

void
MPCPermafrost4::initialize(const Teuchos::Ptr<State>& S) {
  // initialize coupling terms
  S->GetFieldData("surface_subsurface_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_flux", name_)->set_initialized();
  S->GetFieldData("surface_subsurface_energy_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_energy_flux", name_)->set_initialized();

  // Initialize all sub PKs.
  MPCSubsurface::initialize(S);

  // ensure continuity of ICs... surface takes precedence.
  CopySurfaceToSubsurface(*S->GetFieldData("surface-pressure", surf_flow_pk_->name()),
                          S->GetFieldData("pressure", domain_flow_pk_->name()).ptr());
  CopySurfaceToSubsurface(*S->GetFieldData("surface-temperature", surf_energy_pk_->name()),
                          S->GetFieldData("temperature", domain_energy_pk_->name()).ptr());
}

void
MPCPermafrost4::set_states(const Teuchos::RCP<const State>& S,
                           const Teuchos::RCP<State>& S_inter,
                           const Teuchos::RCP<State>& S_next) {
  MPCSubsurface::set_states(S,S_inter,S_next);
  water_->set_states(S,S_inter,S_next);
}


// Compute the non-linear functional g = g(t,u,udot)
void
MPCPermafrost4::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  solution_to_state(*u_new, S_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->Functional(t_old, t_new, u_old->SubVector(2),
                            u_new->SubVector(2), g->SubVector(2));

  // The residual of the surface flow equation provides the mass flux from
  // subsurface to surface.
  Epetra_MultiVector& source = *S_next_->GetFieldData("surface_subsurface_flux",
          name_)->ViewComponent("cell",false);
  source = *g->SubVector(2)->Data()->ViewComponent("cell",false);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->Functional(t_old, t_new, u_old->SubVector(0),
          u_new->SubVector(0), g->SubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(2)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  // Now that mass fluxes are done, do energy.
  // Evaluate the surface energy residual
  surf_energy_pk_->Functional(t_old, t_new, u_old->SubVector(3),
          u_new->SubVector(3), g->SubVector(3));

  // The residual of the surface energy equation provides the diffusive energy
  // flux from subsurface to surface.
  Epetra_MultiVector& esource =
      *S_next_->GetFieldData("surface_subsurface_energy_flux", name_)
      ->ViewComponent("cell",false);
  esource = *g->SubVector(3)->Data()->ViewComponent("cell",false);

  // Evaluate the subsurface energy residual.
  domain_energy_pk_->Functional(t_old, t_new, u_old->SubVector(1),
          u_new->SubVector(1), g->SubVector(1));

  // All energy fluxes have been taken by the subsurface.
  g->SubVector(3)->Data()->ViewComponent("cell",false)->PutScalar(0.);
}

// -- Apply preconditioner
int MPCPermafrost4::ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
        Teuchos::RCP<TreeVector> Pr) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon application:" << std::endl;

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residuals (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  r_ps");
    vnames.push_back("  r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(r->SubVector(2)->Data().ptr());
    vecs.push_back(r->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Residuals (subsurface):" << std::endl;
    vnames[0] = "  r_p";
    vnames[1] = "  r_T";
    vecs[0] = r->SubVector(0)->Data().ptr();
    vecs[1] = r->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  // make a new TreeVector that is just the subsurface values (by pointer).
  // -- note these const casts are necessary to create the new TreeVector, but
  // since the TreeVector COULD be const (it is only used in a single method,
  // in which it is const), const-correctness is not violated here.  The
  // correct solution would be to have a TV constructor that took const
  // subvectors and made a const TV?
  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector());
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(1)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector());
  domain_Pu_tv->PushBack(Pr->SubVector(0));
  domain_Pu_tv->PushBack(Pr->SubVector(1));

  // call the operator's inverse
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying coupled subsurface operator." << std::endl;
  int ierr = linsolve_preconditioner_->ApplyInverse(*domain_u_tv, *domain_Pu_tv);

  // rescale to Pa from MPa
  Pr->SubVector(0)->Data()->Scale(1.e6);

  // Copy subsurface face corrections to surface cell corrections
  CopySubsurfaceToSurface(*Pr->SubVector(0)->Data(),
                          Pr->SubVector(2)->Data().ptr());
  CopySubsurfaceToSurface(*Pr->SubVector(1)->Data(),
                          Pr->SubVector(3)->Data().ptr());

  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "PC * residuals (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC * r_ps");
    vnames.push_back("  PC * r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(Pr->SubVector(2)->Data().ptr());
    vecs.push_back(Pr->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);
  }
  return (ierr > 0) ? 0 : 1;
}

// -- Update the preconditioner.
void
MPCPermafrost4::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update the various components -- note it is important that subsurface are
  // done first (which is handled as they are listed first)
  MPCSubsurface::UpdatePreconditioner(t, up, h, false);
  
  // Add the surface off-diagonal blocks.
  // -- surface dWC/dT
  // -- ALWAYS ZERO
  // S_next_->GetFieldEvaluator("surface-water_content")
  //     ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface-temperature");
  // Teuchos::RCP<const CompositeVector> dWCdT =
  //     S_next_->GetFieldData("dsurface-water_content_dsurface-temperature");
  // dWC_dT_surf_->AddAccumulationTerm(*dWCdT->ViewComponent("cell", false), h);
  
  // -- surface dE_dp
  S_next_->GetFieldEvaluator("surface-energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface-pressure");
  Teuchos::RCP<const CompositeVector> dEdp =
      S_next_->GetFieldData("dsurface-energy_dsurface-pressure");
  dE_dp_surf_->AddAccumulationTerm(*dEdp->ViewComponent("cell", false), h);

  // write for debugging
  std::vector<std::string> vnames;
  //  vnames.push_back("  dwc_dT");
  vnames.push_back("  de_dp");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  //  vecs.push_back(dWCdT.ptr());
  vecs.push_back(dEdp.ptr());
  surf_db_->WriteVectors(vnames, vecs, false);

  // assemble
  // -- scale the pressure dofs
  CompositeVector scaling(sub_pks_[0]->preconditioner()->DomainMap());
  scaling.PutScalar(1.e6); // dWC/dp_Pa * (Pa / MPa) --> dWC/dp_MPa
  sub_pks_[0]->preconditioner()->Rescale(scaling);
  dE_dp_block_->Rescale(scaling);
  
  preconditioner_->AssembleMatrix();

  if (dump_) {
    std::stringstream filename;
    filename << "FullyCoupled_PC_" << S_next_->cycle() << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename.str().c_str(), *preconditioner_->A());
  }

  preconditioner_->InitPreconditioner(plist_->sublist("preconditioner"));
  
  // update ewc Precons if needed
  //  surf_ewc_->UpdatePreconditioner(t, up, h);
}

// -- Modify the predictor.
bool
MPCPermafrost4::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  bool modified = false;

  // write predictor
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Extrapolated Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_extrap");
    vnames.push_back("  Ts_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Extrapolated Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_extrap";
    vnames[1] = "  T_extrap";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  // Make a new TreeVector that is just the subsurface values (by pointer).
  Teuchos::RCP<TreeVector> sub_u = Teuchos::rcp(new TreeVector());
  sub_u->PushBack(u->SubVector(0));
  sub_u->PushBack(u->SubVector(1));

  // Subsurface EWC, modifies cells
  modified |= ewc_->ModifyPredictor(h,sub_u);

  // write predictor
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "EWC Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_extrap");
    vnames.push_back("  Ts_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "EWC Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_extrap";
    vnames[1] = "  T_extrap";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
  
  // Calculate consistent faces
  modified |= domain_flow_pk_->ModifyPredictor(h, u0->SubVector(0), u->SubVector(0));
  modified |= domain_energy_pk_->ModifyPredictor(h, u0->SubVector(1), u->SubVector(1));

  // write predictor
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "EWC/Consistent Face Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_extrap");
    vnames.push_back("  Ts_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "EWC/Consistent Face Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_extrap";
    vnames[1] = "  T_extrap";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
  
  // Copy consistent faces to surface
  if (modified) {
    S_next_->GetFieldEvaluator("surface-relative_permeability")->HasFieldChanged(S_next_.ptr(),name_);
    Teuchos::RCP<const CompositeVector> h_prev = S_inter_->GetFieldData("ponded_depth");
    MergeSubsurfaceAndSurfacePressure(*h_prev, u->SubVector(0)->Data().ptr(), u->SubVector(2)->Data().ptr());
    CopySubsurfaceToSurface(*u->SubVector(1)->Data(), u->SubVector(3)->Data().ptr());
  }

  // Hack surface faces
  bool newly_modified = false;
  newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
  newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
  newly_modified |= water_->ModifyPredictor_TempFromSource(h, u);
  modified |= newly_modified;

  // write predictor
  if (newly_modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Spurt Fixed Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_extrap");
    vnames.push_back("  Ts_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Spurt Fixed Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_extrap";
    vnames[1] = "  T_extrap";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
  

  // -- copy surf --> sub
  //  if (newly_modified) {
  CopySurfaceToSubsurface(*u->SubVector(2)->Data(), u->SubVector(0)->Data().ptr());
  CopySurfaceToSubsurface(*u->SubVector(3)->Data(), u->SubVector(1)->Data().ptr());
  //  }

  // Calculate consistent surface faces
  surf_flow_pk_->ChangedSolution();
  surf_energy_pk_->ChangedSolution();
  modified |= surf_flow_pk_->ModifyPredictor(h, u0->SubVector(2), u->SubVector(2));
  modified |= surf_energy_pk_->ModifyPredictor(h, u0->SubVector(3), u->SubVector(3));

  return modified;
}

// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCPermafrost4::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> r,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();
  
  // dump NKAd correction to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "NKA * PC * residuals (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  NKA*PC*r_ps");
    vnames.push_back("  NKA*PC*r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(du->SubVector(2)->Data().ptr());
    vecs.push_back(du->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "NKA * PC * residuals (subsurface):" << std::endl;
    vnames[0] = "  NKA*PC*r_p";
    vnames[1] = "  NKA*PC*r_T";
    vecs[0] = du->SubVector(0)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  int n_modified = 0;

  // apply limiters
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult base_result =   
    MPCSubsurface::ModifyCorrection(h,r,u,du);
  if (base_result) n_modified++;
  
  // modify correction using water approaches
  n_modified += water_->ModifyCorrection_WaterFaceLimiter(h, r, u, du);
  double damping = water_->ModifyCorrection_WaterSpurtDamp(h, r, u, du);
  n_modified += water_->ModifyCorrection_WaterSpurtCap(h, r, u, du, damping);

  // -- accumulate globally
  int n_modified_l = n_modified;
  u->SubVector(0)->Data()->Comm().SumAll(&n_modified_l, &n_modified, 1);
  bool modified = (n_modified > 0) || (damping < 1.);

  // -- calculate consistent subsurface cells
  // if (modified) {
  //   if (consistent_cells_) {
  //     // Derive subsurface cell corrections.
  //     pc_flow_->UpdateConsistentCellCorrection(
  //         *u->SubVector(0)->Data(),
  //         du->SubVector(0)->Data().ptr());
  //   }
  // }

  // // modify correction on subsurface cells using EWC
  // if (precon_type_ == PRECON_EWC) {
  //     // make a new TreeVector that is just the subsurface values (by pointer).
  //     // -- note these const casts are necessary to create the new TreeVector, but
  //     //    since the TreeVector COULD be const (it is only used in a single method,
  //     //    in which it is const), const-correctness is not violated here.
  //     Teuchos::RCP<TreeVector> domain_res_tv = Teuchos::rcp(new TreeVector());
  //     domain_res_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(0)));
  //     domain_res_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(1)));
      
  //     Teuchos::RCP<TreeVector> domain_du_tv = Teuchos::rcp(new TreeVector());
  //     domain_du_tv->PushBack(du->SubVector(0));
  //     domain_du_tv->PushBack(du->SubVector(1));

  //     // make sure we can back-calc face corrections that preserve residuals on faces
  //     Teuchos::RCP<TreeVector> res0 = Teuchos::rcp(new TreeVector(*domain_res_tv));
  //     res0->PutScalar(0.);
  //     Teuchos::RCP<TreeVector> du_std = Teuchos::rcp(new TreeVector(*domain_du_tv));
  //     *du_std = *domain_du_tv;

  //     // call EWC, which does Pu_p <-- Pu_p_std + dPu_p
  //     sub_ewc_->ApplyPreconditioner(domain_res_tv, domain_du_tv);

  //     // calculate dPu_lambda from dPu_p
  //     du_std->Update(1.0, *domain_du_tv, -1.0);
  //     preconditioner_->UpdateConsistentFaceCorrection(*res0, du_std.ptr());

  //     // update Pu_lambda <-- Pu_lambda_std + dPu_lambda
  //     domain_du_tv->SubVector(0)->Data()->ViewComponent("face",false)->Update(1.,
  //             *du_std->SubVector(0)->Data()->ViewComponent("face",false), 1.);
  //     domain_du_tv->SubVector(1)->Data()->ViewComponent("face",false)->Update(1.,
  //             *du_std->SubVector(1)->Data()->ViewComponent("face",false), 1.);
  //     modified = true;
  // }
  
  if (modified) {
    // Copy subsurface face corrections to surface cell corrections
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
                            du->SubVector(2)->Data().ptr());

    // Derive surface face corrections.
    //    UpdateConsistentFaceCorrectionWater_(r.ptr(), u.ptr(), du.ptr());
  }

  // dump modified correction to screen
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Modified correction:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  Mod NKA*PC*r_ps");
    vnames.push_back("  Mod NKA*PC*r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(du->SubVector(2)->Data().ptr());
    vecs.push_back(du->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Modified correction:" << std::endl;
    vnames[0] = "  Mod NKA*PC*r_p";
    vnames[1] = "  Mod NKA*PC*r_T";
    vecs[0] = du->SubVector(0)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  return modified ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING :
      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

// void
// MPCPermafrost4::UpdateConsistentFaceCorrectionWater_(const Teuchos::Ptr<const TreeVector>& r,
//         const Teuchos::Ptr<const TreeVector>& u,
//         const Teuchos::Ptr<TreeVector>& du) {
//   Teuchos::OSTab tab = vo_->getOSTab();

//   // pull out corrections
//   Teuchos::RCP<CompositeVector> surf_dp = du->SubVector(2)->Data();
//   Epetra_MultiVector& surf_dp_c = *surf_dp->ViewComponent("cell",false);
//   Teuchos::RCP<CompositeVector> surf_dT = du->SubVector(3)->Data();
//   Epetra_MultiVector& surf_dT_c = *surf_dT->ViewComponent("cell",false);

//   // Calculate resulting delta h on the surface
//   Teuchos::RCP<CompositeVector> surf_dh = Teuchos::rcp(new CompositeVector(*surf_dp));
//   surf_dh->PutScalar(0.);

//   // Grab the necessary components of the solution, in both CV and TV forms.
//   Teuchos::RCP<CompositeVector> cv_p = S_next_->GetFieldData("surface-pressure", sub_pks_[2]->name());
//   Teuchos::RCP<CompositeVector> cv_T = S_next_->GetFieldData("surface-temperature", sub_pks_[3]->name());

//   Teuchos::RCP<const TreeVector> tv_p;
//   Teuchos::RCP<const TreeVector> tv_T;
//   if (u == Teuchos::null) {
//     // This method is called by the ApplyPreconditioner(), which does not have access to u.
//     Teuchos::RCP<TreeVector> tv_p_tmp = Teuchos::rcp(new TreeVector());
//     tv_p_tmp->SetData(cv_p);
//     tv_p = tv_p_tmp;

//     Teuchos::RCP<TreeVector> tv_T_tmp = Teuchos::rcp(new TreeVector());
//     tv_T_tmp->SetData(cv_T);
//     tv_T = tv_T_tmp;
//   } else {
//     // This method is called by ModifyCorrection(), which does have access to u.
//     tv_p = u->SubVector(2);
//     tv_T = u->SubVector(3);
//   }

//   // Calculate/get old ponded depth, first ensuring PC didn't result in inadmissible solution
//   if (sub_pks_[2]->IsAdmissible(tv_p) &&
//       sub_pks_[3]->IsAdmissible(tv_T)) {
//     S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
//     *surf_dh->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

//     // new ponded depth (again, checking admissibility)
//     cv_p->ViewComponent("cell",false)->Update(-1., surf_dp_c, 1.);
//     cv_T->ViewComponent("cell",false)->Update(-1., surf_dT_c, 1.);
//     sub_pks_[2]->ChangedSolution();
//     sub_pks_[3]->ChangedSolution();

//     if (sub_pks_[2]->IsAdmissible(tv_p) &&
//         sub_pks_[3]->IsAdmissible(tv_T)) {
//       S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);

//       // put delta ponded depth into surf_dh_cell
//       surf_dh->ViewComponent("cell",false)
//           ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

//       // update delta faces
//       pc_surf_flow_->UpdateConsistentFaceCorrection(*r->SubVector(2)->Data(), surf_dh.ptr());
//       *surf_dp->ViewComponent("face",false) = *surf_dh->ViewComponent("face",false);

//       // dump to screen
//       if (vo_->os_OK(Teuchos::VERB_HIGH)) {
//         *vo_->os() << "Precon water correction." << std::endl;

//         std::vector<std::string> vnames;
//         vnames.push_back("pd_new");
//         vnames.push_back("delta_pd");

//         std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
//         vecs.push_back(S_next_->GetFieldData("ponded_depth").ptr());
//         vecs.push_back(surf_dh.ptr());
//         surf_db_->WriteVectors(vnames, vecs, true);
//       }
//     }

//     // revert solution so we don't break things
//     S_next_->GetFieldData("surface-pressure",sub_pks_[2]->name())
//         ->ViewComponent("cell",false)->Update(1., surf_dp_c, 1.);
//     sub_pks_[2]->ChangedSolution();
//     S_next_->GetFieldData("surface-temperature",sub_pks_[3]->name())
//         ->ViewComponent("cell",false)->Update(1., surf_dT_c, 1.);
//     sub_pks_[3]->ChangedSolution();
//   }
// }


int
MPCPermafrost4::ModifyCorrection_FrozenSurface_(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  const Epetra_MultiVector& res_surf_c = *res->SubVector(2)->Data()->ViewComponent("cell",false);
  const Epetra_MultiVector& u_surf_c = *u->SubVector(2)->Data()->ViewComponent("cell",false);

  Epetra_MultiVector& du_surf_c = *du->SubVector(2)->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& du_domain_f = *du->SubVector(0)->Data()->ViewComponent("face",false);

  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  S_next_->GetFieldEvaluator("surface-water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface-pressure");
  const Epetra_MultiVector& dWC_dp_surf_c =
      *S_next_->GetFieldData("dsurface-water_content_dsurface-pressure")
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
