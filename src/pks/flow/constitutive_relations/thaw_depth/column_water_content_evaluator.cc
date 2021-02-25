/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column water content evaluator gets the subsurface water content and map to surface_star 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "column_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {



ColumnWaterContentEvaluator::ColumnWaterContentEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  auto pos = domain_.find_last_of(':');
  int col_id = std::stoi(domain_.substr(pos+1, domain_.size()));
  
  std::stringstream domain_ss;
  domain_ss << "column:"<< col_id;
  wc_key_ = Keys::getKey(domain_ss.str(),"water_content");
  dependencies_.insert(wc_key_);
}
  

ColumnWaterContentEvaluator::ColumnWaterContentEvaluator(const ColumnWaterContentEvaluator& other)
  : SecondaryVariableFieldEvaluator(other),
    wc_key_(other.wc_key_)
{}
  
Teuchos::RCP<FieldEvaluator>
ColumnWaterContentEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnWaterContentEvaluator(*this));
}


void
ColumnWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  // search through the column and find the deepest unfrozen cell

  const auto& wc_c = *S->GetFieldData(wc_key_)->ViewComponent("cell", false);

  int col_cells = wc_c.MyLength();
  double col_sum = 0;
  
  for (int i=0; i!=col_cells; ++i) {
    col_sum += wc_c[0][i];
  }
  
  res_c[0][0] = col_sum;
 
}
  
void
ColumnWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnWaterContentEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
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
ColumnWaterContentEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{

  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
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
