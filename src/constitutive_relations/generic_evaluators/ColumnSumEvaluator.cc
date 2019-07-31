/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes the thaw depth 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "ColumnSumEvaluator.hh"

namespace Amanzi {
namespace Relations {



ColumnSumEvaluator::ColumnSumEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  std::string name;
  if (!plist.isParameter("evaluator dependency")) {
    if (plist.isParameter("evaluator dependency suffix")) {
      name = plist_.get<std::string>("evaluator dependency suffix");
      Key domain = Keys::getDomain(my_key_);
      Key varname = Keys::getKey(domain, name);
      dependencies_.insert(varname);
      Key pname = name + std::string(" coefficient");
      coefs_[varname] = plist.get<double>(pname, 1.0);
    } else {
      Errors::Message msg;
      msg << "ColumnSumEvaluator for: \"" << my_key_ << "\" has no dependencies.";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    name = plist.get<std::string>("evaluator dependency");
    //    AMANZI_ASSERT( )
    Key pname = name + std::string(" coefficient");
    coefs_[name] = plist.get<double>(pname, 1.0);
  }

  dep_key_ = plist.get<std::string>("evaluator dependency");
  Key domain_name = Keys::getDomain(dep_key_);
  Key surf_name = Keys::getDomain(my_key_);
  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
  surf_cv_key_ = Keys::readKey(plist_, surf_name, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);
  mdl_key_ = Keys::readKey(plist_, domain_name, "molar density", "molar_density_liquid");
  dependencies_.insert(mdl_key_);
}
  

ColumnSumEvaluator::ColumnSumEvaluator(const ColumnSumEvaluator& other) : 
  SecondaryVariableFieldEvaluator(other),
  coefs_(other.coefs_),dep_key_(other.dep_key_),cv_key_(other.cv_key_),surf_cv_key_(other.surf_cv_key_),mdl_key_(other.mdl_key_) {}
  
Teuchos::RCP<FieldEvaluator>
ColumnSumEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnSumEvaluator(*this));
}


void
ColumnSumEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 

  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& dep_c = *S->GetFieldData(dep_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& mld = *S->GetFieldData(mdl_key_)->ViewComponent("cell",false);

  Key domain = Keys::getDomain(my_key_);
  assert(!domain.empty());

  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh = S->GetMesh();
  for (int c=0; c!=res_c.MyLength(); ++c) {
    double sum = 0;
    for (auto i : subsurf_mesh->cells_of_column(c)) {
      sum += dep_c[0][i]*cv[0][i] / (mld[0][i] * surf_cv[0][c]);
    }
    res_c[0][c] = sum;
  }
 
}


void
ColumnSumEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(coefs_[wrt_key]);
}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnSumEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ColumnSumEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{  
  Key domain = Keys::getDomain(my_key_);
  AMANZI_ASSERT(domain == "surface_star" || domain == "surface");
  
  int ncells = S->GetMesh(domain)->num_entities(AmanziMesh::CELL,
          AmanziMesh::Parallel_type::OWNED);
  
  if (domain == domain) {
    for (int c =0; c < ncells; c++){
      //      std::stringstream name;
      //int id = S->GetMesh(domain)->cell_map(false).GID(c);
      //name << "column_"<< id;
      //      Key dep_key = Keys::getKey(domain,dep_key_);
      dependencies_.insert(dep_key_);
    }
  } 
  
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_key_ != std::string(""));
  
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
