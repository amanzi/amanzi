/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#include <algorithm>
#include "boost/algorithm/string/predicate.hpp"

#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "surface_balance_implicit_subgrid.hh"

namespace Amanzi {
namespace SurfaceBalance {

ImplicitSubgrid::ImplicitSubgrid(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  SurfaceBalanceBase(pk_tree, global_list,  S, solution)
{
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "snow_water_equivalent");

  Teuchos::ParameterList& FElist = S->FEList();

  // set up keys
  snow_dens_key_ = Keys::readKey(*plist_, domain_, "snow density", "density");
  snow_age_key_ = Keys::readKey(*plist_, domain_, "snow age", "age");
  new_snow_key_ = Keys::readKey(*plist_, domain_, "new snow source", "source");
  
  // set up additional primary variables -- this is very hacky, and can become an evaluator in new-state
  // -- snow density
  Teuchos::ParameterList& snow_dens_sublist = FElist.sublist(snow_dens_key_);
  snow_dens_sublist.set("evaluator name", snow_dens_key_);
  snow_dens_sublist.set("field evaluator type", "primary variable");
  S->FEList().set(snow_dens_key_, snow_dens_sublist);

  // set the error tolerance for snow
  plist_->set("absolute error tolerance", 0.01);
}

// main methods
// -- Setup data.
void
ImplicitSubgrid::Setup(const Teuchos::Ptr<State>& S) {
  SurfaceBalanceBase::Setup(S);

  // requireiments: things I use
  S->RequireField(new_snow_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(new_snow_key_);
  
  // requirements: other primary variables
  Teuchos::RCP<FieldEvaluator> fm;
  S->RequireField(snow_dens_key_, name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(snow_dens_key_);
  fm = S->GetFieldEvaluator(snow_dens_key_);
  pvfe_snow_dens_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_snow_dens_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  // requirements: internally we must track snow age, this should become an
  // evaluator with snow_density.
  S->RequireField(snow_age_key_, name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
}

// -- Initialize owned (dependent) variables.
void
ImplicitSubgrid::Initialize(const Teuchos::Ptr<State>& S) {
  SurfaceBalanceBase::Initialize(S);

  // initialize snow density, age
  AMANZI_ASSERT(plist_->isSublist("initial condition"));
  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");
  
  if (!S->GetField(snow_dens_key_)->initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S->GetField(snow_dens_key_, name_)->Initialize(ic_list);
      S->GetField(snow_dens_key_, name_)->set_initialized();
      S->GetField(snow_age_key_, name_)->Initialize(ic_list);
      S->GetField(snow_age_key_, name_)->set_initialized();
    } else {
      // initialize density to fresh powder, age to 0
      SEBPhysics::ModelParams params;
      S->GetFieldData(snow_dens_key_,name_)->PutScalar(params.density_freshsnow);
      S->GetField(snow_dens_key_, name_)->set_initialized();
      S->GetFieldData(snow_age_key_,name_)->PutScalar(0.);
      S->GetField(snow_age_key_, name_)->set_initialized();
    }
  }
}


bool
ImplicitSubgrid::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  Epetra_MultiVector& u_vec = *u->Data()->ViewComponent("cell",false);
  unsigned int ncells = u_vec.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    u_vec[0][c] = std::max(0., u_vec[0][c]);
  }
  return true;
}


// computes the non-linear functional g = g(t,u,udot)
void
ImplicitSubgrid::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
        Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  int cycle = S_next_->cycle();

  AMANZI_ASSERT(S_next_->GetFieldEvaluator(new_snow_key_).get() == S_next_->GetFieldEvaluator(source_key_).get());
  
  // update the residual
  SurfaceBalanceBase::FunctionalResidual(t_old, t_new, u_old, u_new, g);  

  // first fill the role of age/density evaluator, as these depend upon old and new values
  const auto& snow_swe_old = *S_inter_->GetFieldData(conserved_key_)->ViewComponent("cell",false);
  const auto& cell_vol = *S_next_->GetFieldData(cell_vol_key_)->ViewComponent("cell",false);

  const auto& snow_age_old = *S_inter_->GetFieldData(snow_age_key_)->ViewComponent("cell",false);
  auto& snow_age_new = *S_next_->GetFieldData(snow_age_key_, name_)->ViewComponent("cell",false);

  const auto& snow_dens_old = *S_inter_->GetFieldData(snow_dens_key_)->ViewComponent("cell",false);
  auto& snow_dens_new = *S_next_->GetFieldData(snow_dens_key_, name_)->ViewComponent("cell",false);
  
  S_next_->GetFieldEvaluator(new_snow_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const auto& new_snow = *S_next_->GetFieldData(new_snow_key_)->ViewComponent("cell",false);

  S_next_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const auto& source = *S_next_->GetFieldData(source_key_)->ViewComponent("cell",false);

  SEBPhysics::ModelParams params;
  double dt_days = (t_new - t_old) / 86400.;
  for (unsigned int c=0; c!=snow_dens_new.MyLength(); ++c) {
    double swe_added = new_snow[0][c] * (t_new - t_old) * cell_vol[0][c];
    double swe_lost = (new_snow[0][c] - source[0][c] ) * (t_new - t_old) * cell_vol[0][c];
    double swe_old = snow_swe_old[0][c];

    double age_new_snow = dt_days / 2.;
    if (swe_old + swe_added - swe_lost < 1.e-10) {
      snow_age_new[0][c] = 0.;
      snow_dens_new[0][c] = params.density_freshsnow;
    } else {
      // age the old snow
      double age_settled = snow_age_old[0][c] + dt_days;
      double dens_settled = params.density_freshsnow * std::max(std::pow(age_settled, 0.3), 1.);
      
      // Match frost age with assigned density -- Calculate which day frost
      // density matched snow defermation function from (Martinec, 1977)
      //double age_frost = std::pow((params.density_frost / params.density_freshsnow),
      //        (1/0.3)) - 1 + dt;

      // ignoring frost, just weighting precip and old snow
      snow_age_new[0][c] = (age_settled * std::max(swe_old - swe_lost,0.) + age_new_snow * swe_added)
                           / (std::max(swe_old - swe_lost,0.) + swe_added);
      snow_dens_new[0][c] = (dens_settled * std::max(swe_old - swe_lost,0.) + params.density_freshsnow * swe_added)
                            / (std::max(swe_old - swe_lost,0.) + swe_added);
    }
  }
  pvfe_snow_dens_->SetFieldAsChanged(S_next_.ptr());


  // debugging
  std::vector<std::string> vnames;
  vnames.push_back("snow age");
  vnames.push_back("snow dens");

  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_next_->GetFieldData(snow_age_key_).ptr());
  vecs.push_back(S_next_->GetFieldData(snow_dens_key_).ptr());
  db_->WriteVectors(vnames, vecs, false);
  
}


} // namespace
} // namespace
