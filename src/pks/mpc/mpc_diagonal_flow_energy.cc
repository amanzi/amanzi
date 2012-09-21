/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "field_evaluator.hh"
#include "mpc_diagonal_flow_energy.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCDiagonalFlowEnergy> MPCDiagonalFlowEnergy::reg_("energy-flow diagonally coupled");


MPCDiagonalFlowEnergy::MPCDiagonalFlowEnergy(Teuchos::ParameterList& mpc_plist,
        const Teuchos::RCP<State>& S, const Teuchos::RCP<TreeVector>& soln) :
    StrongMPC(mpc_plist, S, soln) {};


void MPCDiagonalFlowEnergy::initialize(const Teuchos::RCP<State>& S) {
        
  StrongMPC::initialize(S);

  std::cout<<"before preconditioning approac\n";
  exit(0);
  std::string methodstring = mpc_plist_.get<std::string>("preconditioning approach");
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

  damping_ = mpc_plist_.get<double>("preconditioner damping", 1.0);

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
