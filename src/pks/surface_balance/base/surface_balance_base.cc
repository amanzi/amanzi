/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#include "surface_balance_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

  SurfaceBalanceBase::SurfaceBalanceBase(Teuchos::Ptr<State> S,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           Teuchos::ParameterList& FElist,
           const Teuchos::RCP<TreeVector>& solution) :
    PKPhysicalBDFBase(S, plist, FElist, solution),
    PKDefaultBase(S, plist, FElist, solution)
{
  // name the layer
  layer_ = plist->get<std::string>("layer name", name_);
  source_key_  = getKey(layer_, "source_sink");
  source_key_ = plist->get<std::string>("source key", source_key_);

  theta_ = plist->get<double>("time discretization theta", 1.0);
  ASSERT(theta_ <= 1.);
  ASSERT(theta_ >= 0.);

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01 * 55000.); // h * nl

}

// main methods
// -- Setup data.
void
SurfaceBalanceBase::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);

  // requirements: source terms from above
  S->RequireField(source_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(source_key_);

  conserved_quantity_ = conserved_key_ != key_;
  if (conserved_quantity_) {
    S->RequireField(conserved_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(conserved_key_);
  }    
}


// computes the non-linear functional g = g(t,u,udot)
void
SurfaceBalanceBase::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                            Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();
  double dt = t_new - t_old;
  double T_eps = 0.0001;

  bool debug = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) debug = true;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl;
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("u_new"); vnames.push_back("u_old");
    vecs.push_back(u_new->Data().ptr());
    vecs.push_back(u_old->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();
  }

  S_next_->GetFieldEvaluator(cell_vol_key_)->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> cv = S_next_->GetFieldData(cell_vol_key_);

  if (conserved_quantity_) {
    S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> conserved1 = S_next_->GetFieldData(conserved_key_);
    Teuchos::RCP<const CompositeVector> conserved0 = S_inter_->GetFieldData(conserved_key_);
    g->Data()->Update(1.0/dt, *conserved1, -1.0/dt, *conserved0, 0.0);
  } else {
    g->Update(1.0/dt, *u_new, -1.0/dt, *u_old, 0.0);
  }
  db_->WriteVector("res(acc)", g->Data().ptr());
  //  db_->WriteVector("  sM", S_next_->GetFieldData("macropore_saturation_liquid").ptr());

  if (theta_ < 1.0) {
    S_inter_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_inter_.ptr(), name_);
    g->Data()->Multiply(-(1.0 - theta_), *S_inter_->GetFieldData(source_key_), *cv, 1.);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      db_->WriteVector("source0", S_inter_->GetFieldData(source_key_).ptr(), false);
      //      db_->WriteVector(" drainage0", S_inter_->GetFieldData("litter_drainage").ptr(), false);
    }
  }
  if (theta_ > 0.0) {
    S_next_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);
    g->Data()->Multiply(-theta_, *S_next_->GetFieldData(source_key_), *cv, 1.);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      db_->WriteVector("source1", S_next_->GetFieldData(source_key_).ptr(), false);
      //      db_->WriteVector(" drainage1", S_next_->GetFieldData("litter_drainage").ptr(), false);
    }
  }
  db_->WriteVector("res(source)", g->Data().ptr());
}

  
// updates the preconditioner
void
SurfaceBalanceBase::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(*up, S_next_);

  if (conserved_quantity_) {
    if (jac_ == Teuchos::null) {
      jac_ = Teuchos::rcp(new CompositeVector(*up->Data()));
    }

    S_next_->GetFieldEvaluator(conserved_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    std::string dkey = std::string("d")+conserved_key_+std::string("_d")+key_;
    *jac_ = *S_next_->GetFieldData(dkey);
    jac_->Scale(1./h);

    if (S_next_->GetFieldEvaluator(source_key_)->IsDependency(S_next_.ptr(), key_)) {
      S_next_->GetFieldEvaluator(source_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
      std::string dkey = std::string("d")+source_key_+std::string("_d")+key_;
      jac_->Multiply(-theta_, *S_next_->GetFieldData(dkey),
                     *S_next_->GetFieldData(cell_vol_key_), 1.);
    }      
  }
}


// applies preconditioner to u and returns the result in Pu
int SurfaceBalanceBase::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

  if (conserved_quantity_) {
    db_->WriteVector("seb_res", u->Data().ptr(), true);
    Pu->Data()->ReciprocalMultiply(1., *jac_, *u->Data(), 0.);
    db_->WriteVector("jac", jac_.ptr(), true);
    db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);
  } else {
    *Pu = *u;
    Pu->Scale(S_next_->time() - S_->time());
  }
  
  return 0;
}


} // namespace
} // namespace
