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

SurfaceBalanceBase::SurfaceBalanceBase(
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           Teuchos::ParameterList& FElist,
           const Teuchos::RCP<TreeVector>& solution) :
    PKPhysicalBDFBase(plist, FElist, solution),
    PKDefaultBase(plist, FElist, solution)
{
  // name the layer
  layer_ = plist->get<std::string>("layer name", name_);
  source_key_  = layer_+std::string("_source");
  source_key_ = plist->get<std::string>("source key", source_key_);

  if (plist->isParameter("conserved quantity")) {
    conserved_quantity_ = true;
    conserved_key_ = plist->get<std::string>("conserved quantity");
  } 

  theta_ = plist->get<double>("time discretization theta", 0.5);
  ASSERT(theta_ <= 1.);
  ASSERT(theta_ >= 0.);
  
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

  S_next_->GetFieldEvaluator("surface_cell_volume")->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> cv = S_next_->GetFieldData("surface_cell_volume");

  if (conserved_quantity_) {
    S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> conserved1 = S_next_->GetFieldData(conserved_key_);
    Teuchos::RCP<const CompositeVector> conserved0 = S_inter_->GetFieldData(conserved_key_);
    g->Data()->Update(1.0/dt, *conserved1, -1.0/dt, *conserved0, 0.0);
  } else {
    g->Update(1.0/dt, *u_new, -1.0/dt, *u_old, 0.0);
  }

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
}

  



} // namespace
} // namespace
