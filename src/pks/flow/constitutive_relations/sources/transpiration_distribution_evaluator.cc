/*
  Distributes transpiration based upon a rooting depth and a wilting-point water-potential factor.
*/

#include "transpiration_distribution_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
TranspirationDistributionEvaluator::TranspirationDistributionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  InitializeFromPlist_();
}


// Copy constructor
TranspirationDistributionEvaluator::TranspirationDistributionEvaluator(const TranspirationDistributionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    p_key_(other.p_key_),
    f_root_key_(other.f_root_key_),
    trans_total_key_(other.trans_total_key_),    
    cv_key_(other.cv_key_),
    surf_cv_key_(other.surf_cv_key_),
    wp_min_(other.wp_min_)
{}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
TranspirationDistributionEvaluator::Clone() const
{
  return Teuchos::rcp(new TranspirationDistributionEvaluator(*this));
}


// Initialize by setting up dependencies
void
TranspirationDistributionEvaluator::InitializeFromPlist_()
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
    Errors::Message message("TranspirationDistributionEvaluator: subsurface domain name and surface domain name cannot be the same.");
    Exceptions::amanzi_throw(message);
  }

  // - pull Keys from plist
  // dependency: pressure
  p_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(p_key_);

  // dependency: rooting_depth_fraction
  f_root_key_ = Keys::readKey(plist_, domain_name, "rooting depth fraction", "rooting_depth_fraction");
  dependencies_.insert(f_root_key_);

  // dependency: transpiration
  trans_total_key_ = Keys::readKey(plist_, surface_name, "total transpiration", "transpiration");
  dependencies_.insert(trans_total_key_);

  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
  surf_cv_key_ = Keys::readKey(plist_, surface_name, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);

  // wilting point
  wp_min_ = plist_.get<double>("wilting point [Pa]", -2.0e6);
}


void
TranspirationDistributionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  // on the subsurface
  const Epetra_MultiVector& p = *S->GetFieldData(p_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& f_root = *S->GetFieldData(f_root_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);

  // on the surface
  const Epetra_MultiVector& trans_total = *S->GetFieldData(trans_total_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);

  double p_atm = *S->GetScalarData("atmospheric_pressure");

  // result, on the subsurface
  const AmanziMesh::Mesh& subsurf_mesh = *result->Mesh();
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);
  
  for (int sc=0; sc!=trans_total.MyLength(); ++sc) {
    double column_total = 0.;
    for (auto c : subsurf_mesh.cells_of_column(sc)) {
      double water_potential_factor = (wp_min_ - p[0][c]) / (wp_min_ - p_atm);
      water_potential_factor = water_potential_factor > 1.0 ? 1 :
                               (water_potential_factor < 0. ? 0. : water_potential_factor);

      result_v[0][c] = trans_total[0][sc] * water_potential_factor * f_root[0][c];
      column_total += result_v[0][c];
    }
    if (column_total <= 0. && trans_total[0][sc] > 0.) {
      Errors::Message message("TranspirationDistributionEvaluator: Broken run, non-zero transpiration draw but no cells with some roots are above the wilting point.");
      Exceptions::amanzi_throw(message);
    }

    if (column_total > 0.) {
      double rescaling_factor = column_total / trans_total[0][sc];
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[0][c] *= rescaling_factor;
      }
    }
  }
}


void
TranspirationDistributionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(0.); // this would be a nontrivial calculation, as it is technically nonlocal due to rescaling issues?
}


void
TranspirationDistributionEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(!my_key_.empty());
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  Key domain = Keys::getDomain(my_key_);
  my_fac->SetMesh(S->GetMesh(domain))
      ->SetComponent("cell", AmanziMesh::CELL, 1)
      ->SetGhosted(true);

  // Create an unowned factory to check my dependencies.
  // -- first those on the subsurface mesh
  Teuchos::RCP<CompositeVectorSpace> dep_fac =
      Teuchos::rcp(new CompositeVectorSpace(*my_fac));
  dep_fac->SetOwned(false);
  S->RequireField(p_key_)->Update(*dep_fac);
  S->RequireField(cv_key_)->Update(*dep_fac);
  S->RequireField(f_root_key_)->Update(*dep_fac);

  // -- next those on the surface mesh
  Teuchos::RCP<CompositeVectorSpace> surf_fac =
      Teuchos::rcp(new CompositeVectorSpace());
  surf_fac->SetMesh(S->GetMesh(Keys::getDomain(surf_cv_key_)))
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(trans_total_key_)->Update(*surf_fac);
  S->RequireField(surf_cv_key_)->Update(*surf_fac);

  // Recurse into the tree to propagate info to leaves.
  for (KeySet::const_iterator key=dependencies_.begin();
       key!=dependencies_.end(); ++key) {
    S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
  }
}


} //namespace
} //namespace
} //namespace
