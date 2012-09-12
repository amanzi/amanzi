/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "field_evaluator.hh"
#include "wrm_richards_evaluator.hh"
#include "mpc_diagonal_flow_energy.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCDiagonalFlowEnergy> MPCDiagonalFlowEnergy::reg_("energy-flow diagonally coupled");


bool MPCDiagonalFlowEnergy::advance(double dt) {
  n_iter_ = 0;
  StrongMPC::advance(dt);
}

void MPCDiagonalFlowEnergy::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
        Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  n_iter_++;
  StrongMPC::fun(t_old, t_new, u_old, u_new, g);
}

void MPCDiagonalFlowEnergy::initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC::initialize(S);

  richards_pk_ = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[0]);
  ASSERT(richards_pk_ != Teuchos::null);
  two_phase_pk_ = Teuchos::rcp_dynamic_cast<Energy::TwoPhase>(sub_pks_[1]);
  ASSERT(two_phase_pk_ != Teuchos::null);
  // this is just to check that we have access to protected/private data in
  // Richards/TwoPhase
  ASSERT(richards_pk_->preconditioner_ != Teuchos::null);
  ASSERT(two_phase_pk_->preconditioner_ != Teuchos::null);

  std::string methodstring = plist_.get<std::string>("preconditioning approach");
  if (methodstring == "diagonal") {
    method_ = PRECON_DIAGONAL;
  } else if (methodstring == "upper triangular") {
    method_ = PRECON_UPPER_TRIANGULAR;
  } else if (methodstring == "lower triangular") {
    method_ = PRECON_LOWER_TRIANGULAR;
  } else if (methodstring == "alternating") {
    method_ = PRECON_ALTERNATING;
  } else if (methodstring == "accumulation") {
    method_ = PRECON_ACCUMULATION;
  } else {
    ASSERT(0);
  }

  damping_ = plist_.get<double>("preconditioner damping", 1.0);

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");

  ASSERT(temp->size("cell") == pres->size("cell"));

  if (method_ == PRECON_UPPER_TRIANGULAR || method_ == PRECON_ALTERNATING) {
    D_pT_ = Teuchos::rcp(new Epetra_MultiVector(*pres->ViewComponent("cell", false)));
  }
  if (method_ == PRECON_LOWER_TRIANGULAR || method_ == PRECON_ALTERNATING) {
    D_Tp_ = Teuchos::rcp(new Epetra_MultiVector(*temp->ViewComponent("cell", false)));
  }
};


bool MPCDiagonalFlowEnergy::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  // modification of the initial guess occurs by updating T using this guess,
  // then calculating the p that would be required to keep the same water
  // mass.

  // Stuff temperature into state
  Teuchos::RCP<TreeVector> uT = u->SubVector("energy");
  sub_pks_[1]->solution_to_state(uT, S_next_);

#if DEBUG_FLAG
  std::cout << std::endl;
  std::cout << "Modifying guess:" << std::endl;
  std::cout << "  T0: " << (*uT->data())("cell",0) << std::endl;
  std::cout << "  T1: " << (*uT->data())("cell",99) << std::endl;
#endif

  // update water content, which will get all the needed vals updated at the new temp.
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldChanged(S_next_.ptr(), "richards_pk");

  // get the needed vals
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");
  Teuchos::RCP<const CompositeVector> cv = S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> phi = S_next_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> n_g = S_next_->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> omega_g = S_next_->GetFieldData("mol_frac_gas");
  Teuchos::RCP<const CompositeVector> n_l = S_next_->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> n_i = S_next_->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> one_on_A = S_next_->GetFieldData("wrm_permafrost_one_on_A");
  Teuchos::RCP<const double> p_atm = S_next_->GetScalarData("atmospheric_pressure");

  // get the WRMs
  Teuchos::RCP<FieldEvaluator> wrm_B_eval_base = S_next_->GetFieldEvaluator("wrm_permafrost_one_on_B");
  Teuchos::RCP<Flow::FlowRelations::WRMRichardsEvaluator> wrm_B_eval =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::WRMRichardsEvaluator>(wrm_B_eval_base);
  Teuchos::RCP<Flow::FlowRelations::WRMRegionPairList> wrms = wrm_B_eval->get_WRMs();


  // get the result pressure
  Teuchos::RCP<CompositeVector> pres = u->SubVector("flow")->data();
  for (Flow::FlowRelations::WRMRegionPairList::iterator region=wrms->begin();
       region!=wrms->end(); ++region) {
    std::string name = region->first;
    int ncells = pres->mesh()->get_set_size(name, AmanziMesh::CELL, AmanziMesh::OWNED);
    std::vector<int> cells(ncells);
    pres->mesh()->get_set_entities(name, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      if ((*pres)("cell",*c) < *p_atm) {
        double A_minus_one = (1.0/(*one_on_A)("cell",*c) - 1.0);
        if (*c==0) std::cout << "   A-1(0) = " << A_minus_one << std::endl;
        if (*c==99) std::cout << "   A-1(99) = " << A_minus_one << std::endl;

        double wc = (*wc0)("cell",*c) / (*cv)("cell",*c);
        double sstar = (wc - (*n_g)("cell",*c)*(*omega_g)("cell",*c)*(*phi)("cell",*c)) /
            ((*phi)("cell",*c) * ((*n_l)("cell",*c) - (*n_g)("cell",*c)*(*omega_g)("cell",*c)
                    + (*n_i)("cell",*c)*A_minus_one) - A_minus_one*wc);

        if (*c==0) std::cout << "   S*(0) = " << sstar << std::endl;
        if (*c==99) std::cout << "   S*(99) = " << sstar << std::endl;

        double pc = region->second->capillaryPressure(sstar);

        if (*c==0) std::cout << "   pc(0) = " << pc << std::endl;
        if (*c==99) std::cout << "   pc(99) = " << pc << std::endl;

        (*pres)("cell",*c) = *p_atm - pc;

        if (*c==0) std::cout << "   p(0) = " << (*pres)("cell",*c) << std::endl;
        if (*c==99) std::cout << "   p(99) = " <<  (*pres)("cell",*c) << std::endl;

      }
    }
  }

  return true;
};


void MPCDiagonalFlowEnergy::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  StrongMPC::update_precon(t, up, h);

  if (method_ == PRECON_UPPER_TRIANGULAR || method_ == PRECON_ALTERNATING) {
    // d pressure_residual / d T
    // -- update the accumulation derivative
    S_next_->GetFieldEvaluator("water_content")
        ->HasFieldDerivativeChanged(S_next_.ptr(), "richards_pk", "temperature");

    // -- get the accumulation deriv
    Teuchos::RCP<const CompositeVector> dwc_dT =
        S_next_->GetFieldData("dwater_content_dtemperature");
    for (int c=0; c!=dwc_dT->size("cell"); ++c) {
      (*D_pT_)[0][c] = (*dwc_dT)("cell",c) / h;
    }
  }

  if (method_ == PRECON_LOWER_TRIANGULAR || method_ == PRECON_ALTERNATING) {
    // d temperature_residual / d p
    // -- update the accumulation derivative
    S_next_->GetFieldEvaluator("energy")
        ->HasFieldDerivativeChanged(S_next_.ptr(), "energy_pk", "pressure");

    // -- get the accumulation deriv
    Teuchos::RCP<const CompositeVector> de_dp =
        S_next_->GetFieldData("denergy_dpressure");
    for (int c=0; c!=de_dp->size("cell"); ++c) {
      (*D_Tp_)[0][c] = (*de_dp)("cell",c) / h;
    }
  }
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void MPCDiagonalFlowEnergy::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  if (method_ == PRECON_DIAGONAL) {
    precon_diagonal(u, Pu);
  } else if (method_ == PRECON_UPPER_TRIANGULAR) {
    precon_upper_triangular(u, Pu);
  } else if (method_ == PRECON_LOWER_TRIANGULAR) {
    precon_lower_triangular(u, Pu);
  } else if (method_ == PRECON_ALTERNATING) {
    precon_alternating(u, Pu);
  } else if (method_ == PRECON_ACCUMULATION) {
    precon_accumulation(u, Pu);
  }

  if (damping_ < 1.0) {
    Pu->Scale(damping_);
  }
};


void MPCDiagonalFlowEnergy::precon_diagonal(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  StrongMPC::precon(u, Pu);
};


void MPCDiagonalFlowEnergy::precon_upper_triangular(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

  // pull out the temperature sub-vector
  Teuchos::RCP<const TreeVector> uT = u->SubVector("energy");
  Teuchos::RCP<TreeVector> PuT = Pu->SubVector("energy");

  // apply temperature inverse
  sub_pks_[1]->precon(uT, PuT);

  // subtract off the d pressure residual dT part.
  Teuchos::RCP<TreeVector> up_minusT = Teuchos::rcp(new TreeVector(*u->SubVector("flow")));
  *up_minusT = *u->SubVector("flow");
  up_minusT->data()->ViewComponent("cell",false)->Multiply(-1.0, *PuT->data()->ViewComponent("cell", false), *D_pT_, 1.0);

  // Apply the pressure inverse
  Teuchos::RCP<TreeVector> Pup = Pu->SubVector("flow");
  sub_pks_[0]->precon(up_minusT, Pup);
};


void MPCDiagonalFlowEnergy::precon_lower_triangular(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // pull out the pressure sub-vector
  Teuchos::RCP<const TreeVector> up = u->SubVector("flow");
  Teuchos::RCP<TreeVector> Pup = Pu->SubVector("flow");

  // apply pressure inverse
  sub_pks_[0]->precon(up, Pup);

  // subtract off the d temp residual dp part.
  Teuchos::RCP<TreeVector> uT_minusp = Teuchos::rcp(new TreeVector(*u->SubVector("energy")));
  *uT_minusp = *u->SubVector("energy");
  uT_minusp->data()->ViewComponent("cell",false)->Multiply(-1.0, *Pup->data()->ViewComponent("cell", false), *D_Tp_, 1.0);

  // Apply the temp inverse
  Teuchos::RCP<TreeVector> PuT = Pu->SubVector("energy");
  sub_pks_[1]->precon(uT_minusp, PuT);
};


void MPCDiagonalFlowEnergy::precon_alternating(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  if (n_iter_%2 == 1) {
    precon_upper_triangular(u,Pu);
  } else {
    precon_lower_triangular(u,Pu);
  }
};


void MPCDiagonalFlowEnergy::precon_accumulation(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

  // Trying some stupid strategies.
  //   -- In the first few iterations, just deal with local adjustments of T/p.
  //   -- In later iterations, do the full PC.
  // get the accumulation derivatives
  Teuchos::RCP<const CompositeVector> de_dp =
      S_next_->GetFieldData("denergy_dpressure");
  Teuchos::RCP<const CompositeVector> de_dT =
      S_next_->GetFieldData("denergy_dtemperature");
  Teuchos::RCP<const CompositeVector> dwc_dp =
      S_next_->GetFieldData("dwater_content_dpressure");
  Teuchos::RCP<const CompositeVector> dwc_dT =
      S_next_->GetFieldData("dwater_content_dtemperature");

  // get h
  double h = S_next_->time() - S_inter_->time();

  // get composite vectors
  Teuchos::RCP<const CompositeVector> u_p = u->SubVector("flow")->data();
  Teuchos::RCP<const CompositeVector> u_T = u->SubVector("energy")->data();
  Teuchos::RCP<CompositeVector> Pu_p = Pu->SubVector("flow")->data();
  Teuchos::RCP<CompositeVector> Pu_T = Pu->SubVector("energy")->data();

  for (int c=0; c!=de_dp->size("cell",false); ++c) {
    // local solves
    double h_on_det = h / ((*dwc_dp)("cell",c)*(*de_dT)("cell",c)
                           - (*dwc_dT)("cell",c)*(*de_dp)("cell",c));
    (*Pu_p)("cell",c) = h_on_det * ( (*de_dT)("cell",c)*(*u_p)("cell",c)
            - (*dwc_dT)("cell",c)*(*u_T)("cell",c));
    (*Pu_T)("cell",c) = h_on_det * ( -(*de_dp)("cell",c)*(*u_p)("cell",c)
            + (*dwc_dp)("cell",c)*(*u_T)("cell",c));
  }

#if DEBUG_FLAG
  std::cout << "Precon application:" << std::endl;
  std::cout << "  T0: " << (*u_T)("cell",0,0) << " " << (*u_T)("face",0,3) << std::endl;
  std::cout << "  T1: " << (*u_T)("cell",0,99) << " " << (*u_T)("face",0,497) << std::endl;
  std::cout << "  PC*T0: " << (*Pu_T)("cell",0,0) << " " << (*Pu_T)("face",0,3) << std::endl;
  std::cout << "  PC*T1: " << (*Pu_T)("cell",0,99) << " " << (*Pu_T)("face",0,497) << std::endl;
  std::cout << "Precon application:" << std::endl;
  std::cout << "  p0: " << (*u_p)("cell",0,0) << " " << (*u_p)("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u_p)("cell",0,99) << " " << (*u_p)("face",0,497) << std::endl;
  std::cout << "  PC*p0: " << (*Pu_p)("cell",0,0) << " " << (*Pu_p)("face",0,3) << std::endl;
  std::cout << "  PC*p1: " << (*Pu_p)("cell",0,99) << " " << (*Pu_p)("face",0,497) << std::endl;
#endif
};

} //namespace
