/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes the thaw depth 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "thaw_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {



ThawDepthEvaluator::ThawDepthEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  std::string domain_name=Keys::getDomain(my_key_);
  my_key_ = plist_.get<std::string>("thaw depth key", Keys::getKey(domain_name, "thaw_depth"));
}
  

ThawDepthEvaluator::ThawDepthEvaluator(const ThawDepthEvaluator& other)
    : SecondaryVariableFieldEvaluator(other)
{}
  
Teuchos::RCP<FieldEvaluator>
ThawDepthEvaluator::Clone() const
{
  return Teuchos::rcp(new ThawDepthEvaluator(*this));
}


void
ThawDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; c++){
    int id = S->GetMesh("surface_star")->cell_map(false).GID(c);
    std::stringstream domain;
    domain << "column_" << id;

    const auto& top_z_centroid = S->GetMesh(domain.str())->face_centroid(0);
    AmanziGeometry::Point z_centroid(top_z_centroid);

    // search through the column and find the deepest unfrozen cell
    const auto& temp_c = *S->GetFieldData(Keys::getKey(domain.str(),"temperature"))
                         ->ViewComponent("cell", false);
    int col_cells = temp_c.MyLength();
    for (int i=0; i!=col_cells; ++i) {
      if (temp_c[0][i] >= 273.25) { // this hard codes in the transition width to 0.2 K
        z_centroid = S->GetMesh(domain.str())->face_centroid(i+1);
      }
    }
    
    res_c[0][c] = top_z_centroid[2] - z_centroid[2];
  }
}
  
void
ThawDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
ThawDepthEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
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
ThawDepthEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{  
  Key domain = Keys::getDomain(my_key_);
  AMANZI_ASSERT(domain == "surface_star");
  
  int ncells = S->GetMesh("surface_star")->num_entities(AmanziMesh::CELL,
          AmanziMesh::Parallel_type::OWNED);
  
  if (domain == "surface_star") {
    for (int c =0; c < ncells; c++){
      std::stringstream name;
      int id = S->GetMesh("surface_star")->cell_map(false).GID(c);
      name << "column_"<< id;
      Key temp_key = Keys::getKey(name.str(),"temperature");
      dependencies_.insert(temp_key);
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
