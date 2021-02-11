/*
  The rooting depth fraction evaluator is an algebraic evaluator of a given model.
Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

  
  Generated via evaluator_generator.
*/

#include "rooting_depth_fraction_evaluator.hh"
#include "rooting_depth_fraction_model.hh"

namespace Amanzi {
namespace LandCover {
namespace Relations {

// Constructor from ParameterList
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  if (!plist_.isSublist("rooting_depth_fraction parameters")) {
    Errors::Message message("RootingDepthFractionEvaluator: changed spec -- now must be list of models to match list of PFTs.");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList& sublist = plist_.sublist("rooting_depth_fraction parameters");
  for (auto p : sublist) {
    if (!sublist.isSublist(p.first)) {
      Errors::Message message("RootingDepthFractionEvaluator: changed spec -- now must be list of models to match list of PFTs.");
      Exceptions::amanzi_throw(message);
    }
    models_.push_back(Teuchos::rcp(new RootingDepthFractionModel(sublist.sublist(p.first))));
  }
  InitializeFromPlist_();
}


// Copy constructor
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(const RootingDepthFractionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    z_key_(other.z_key_),    
    cv_key_(other.cv_key_),    
    surf_cv_key_(other.surf_cv_key_),    
    models_(other.models_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
RootingDepthFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new RootingDepthFractionEvaluator(*this));
}


// Initialize by setting up dependencies
void
RootingDepthFractionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomainPrefix(my_key_);
  Key surface_name;
  if (domain_name.empty() || domain_name == "domain") {
    surface_name = "surface";
  } else {
    surface_name = Key("surface_")+domain_name;
  }
  surface_name = plist_.get<std::string>("surface domain name", surface_name);
  if (surface_name == domain_name) {
    Errors::Message message("RootingDepthFractionEvaluator: subsurface domain name and surface domain name cannot be the same.");
    Exceptions::amanzi_throw(message);
  }

  // - pull Keys from plist
  // dependency: depth
  z_key_ = Keys::readKey(plist_, domain_name, "depth", "depth");
  dependencies_.insert(z_key_);

  // cell volume, surface area
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
  surf_cv_key_ = Keys::readKey(plist_, surface_name, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);

  npfts_ = models_.size();
}


void
RootingDepthFractionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& z = *S->GetFieldData(z_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);
  auto& subsurf_mesh = *result->Mesh();

  for (int pft=0; pft!=models_.size(); ++pft) {
    for (int sc=0; sc!=surf_cv.MyLength(); ++sc) {
      double column_total = 0.;
      double f_root_total = 0.;
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[pft][c] = models_[pft]->RootingDepthFraction(z[0][c]);
        column_total += result_v[pft][c] * cv[0][c];
      }
      
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[pft][c] = result_v[pft][c] * 1.0 * surf_cv[0][sc] / column_total;  // normalize to 1
      }
    }
  }
}


void
RootingDepthFractionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  // this should only change if the mesh deforms.  don't do that!
  result->PutScalar(0.);
  AMANZI_ASSERT(0);
}


void
RootingDepthFractionEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(!my_key_.empty());
  auto my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  Key domain = Keys::getDomain(my_key_);
  my_fac->SetMesh(S->GetMesh(domain))
      ->SetComponent("cell", AmanziMesh::CELL, npfts_)
      ->SetGhosted(true);

  // Create an unowned factory to check my dependencies.
  CompositeVectorSpace dep_fac_one;
  dep_fac_one.SetMesh(my_fac->Mesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(z_key_)->Update(dep_fac_one);
  S->RequireField(cv_key_)->Update(dep_fac_one);

  CompositeVectorSpace surf_fac_one;
  surf_fac_one.SetMesh(S->GetMesh(Keys::getDomain(surf_cv_key_)))
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
