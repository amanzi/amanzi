/*
  The rooting depth fraction evaluator is an algebraic evaluator of a given model.
Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

  
  Generated via evaluator_generator.
*/

#include "onepft_rooting_depth_fraction_evaluator.hh"
#include "rooting_depth_fraction_model.hh"

namespace Amanzi {
namespace LandCover {
namespace Relations {

// Constructor from ParameterList
OnePFTRootingDepthFractionEvaluator::OnePFTRootingDepthFractionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("rooting_depth_fraction parameters");
  for (auto p : sublist) {
    if (!sublist.isSublist(p.first)) {
      Errors::Message message("OnePFTRootingDepthFractionEvaluator: expected list of models.");
      Exceptions::amanzi_throw(message);
    }
    auto& model_plist = sublist.sublist(p.first);
    models_.emplace_back(std::make_pair(model_plist.get<std::string>("region"),
            Teuchos::rcp(new RootingDepthFractionModel(model_plist))));
  }
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
OnePFTRootingDepthFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new OnePFTRootingDepthFractionEvaluator(*this));
}


// Initialize by setting up dependencies
void
OnePFTRootingDepthFractionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  subsurf_domain_ = Keys::getDomainPrefix(my_key_);
  if (subsurf_domain_.empty() || subsurf_domain_ == "domain") {
    surf_domain_ = "surface";
  } else {
    surf_domain_ = Key("surface_")+subsurf_domain_;
  }
  surf_domain_ = plist_.get<std::string>("surface domain name", surf_domain_);
  if (surf_domain_ == subsurf_domain_) {
    Errors::Message message("OnePFTRootingDepthFractionEvaluator: subsurface domain name and surface domain name cannot be the same.");
    Exceptions::amanzi_throw(message);
  }

  // - pull Keys from plist
  // dependency: depth
  z_key_ = Keys::readKey(plist_, subsurf_domain_, "depth", "depth");
  dependencies_.insert(z_key_);

  // cell volume, surface area
  cv_key_ = Keys::readKey(plist_, subsurf_domain_, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
  surf_cv_key_ = Keys::readKey(plist_, surf_domain_, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);
}


void
OnePFTRootingDepthFractionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& z = *S->GetFieldData(z_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);
  auto& subsurf_mesh = *S->GetMesh(subsurf_domain_);
  auto& surf_mesh = *S->GetMesh(surf_domain_);

  for (auto region_model : models_) {
    AmanziMesh::Entity_ID_List surf_cells;
    surf_mesh.get_set_entities(region_model.first, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &surf_cells);
    for (auto sc : surf_cells) {
      double column_total = 0.;
      double f_root_total = 0.;
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[0][c] = region_model.second->RootingDepthFraction(z[0][c]);
        column_total += result_v[0][c] * cv[0][c];
      }
      
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[0][c] = result_v[0][c] * 1.0 * surf_cv[0][sc] / column_total;  // normalize to 1
      }
    }
  }
}


void
OnePFTRootingDepthFractionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  // this should only change if the mesh deforms.  don't do that!
  result->PutScalar(0.);
  AMANZI_ASSERT(0);
}


void
OnePFTRootingDepthFractionEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(!my_key_.empty());
  auto my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  my_fac->SetMesh(S->GetMesh(subsurf_domain_))
      ->SetComponent("cell", AmanziMesh::CELL, 1)
      ->SetGhosted(true);

  // Create an unowned factory to check my dependencies.
  CompositeVectorSpace dep_fac_one;
  dep_fac_one.SetMesh(my_fac->Mesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(z_key_)->Update(dep_fac_one);
  S->RequireField(cv_key_)->Update(dep_fac_one);

  CompositeVectorSpace surf_fac_one;
  surf_fac_one.SetMesh(S->GetMesh(surf_domain_))
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(surf_cv_key_)->Update(surf_fac_one);
  
  // Recurse into the tree to propagate info to leaves.
  for (KeySet::const_iterator key=dependencies_.begin();
       key!=dependencies_.end(); ++key) {
    S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
  }
}


} //namespace
} //namespace
} //namespace
