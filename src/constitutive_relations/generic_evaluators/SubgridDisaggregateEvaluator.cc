/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

 * `"source domain name`" ``[string]`` Domain name of the source mesh.

ONE OF:
* `"field suffix`" ``[string]`` Set the suffix of the variable
OR
* `"field key`" ``[string]`` **DOMAIN-FIELD_SUFFIX** 

  
 */

#include "SubgridDisaggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridDisaggregateEvaluator::SubgridDisaggregateEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    source_gid_(-1)
{
  // my_key_ = "surface_column_6-del_max"
  domain_ = Keys::getDomain(my_key_);  // "surface_column_6"
  source_domain_ = plist_.get<std::string>("source domain name"); // "surface_star"
  source_key_ = Keys::readKey(plist_, source_domain_, "field", Keys::getVarName(my_key_)); // "surface_star-del_max"
  dependencies_.insert(source_key_);
}

Teuchos::RCP<FieldEvaluator>
SubgridDisaggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridDisaggregateEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
SubgridDisaggregateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  if (source_gid_ < 0) {
    auto pos = domain_.find_last_of(':');
    AMANZI_ASSERT(pos != domain_.size());
    AMANZI_ASSERT(pos != domain_.size()-1);
    int col_id = std::stoi(domain_.substr(pos+1, domain_.size()));
    AMANZI_ASSERT(col_id >= 0 && col_id <= 1e3);
    source_gid_ = col_id; // 6
  }

  const auto& source = *S->GetFieldData(source_key_)->ViewComponent("cell",false);

  int source_lid = S->GetMesh("surface")->cell_map(false).LID(source_gid_);
  AMANZI_ASSERT(source.MyLength() > source_lid);
  (*result->ViewComponent("cell", false))[0][0] = source[0][source_lid];
}

void
SubgridDisaggregateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(1.);
}


void
SubgridDisaggregateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  S->RequireField(my_key_, my_key_)
      ->SetMesh(S->GetMesh(domain_))
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // Check plist for vis or checkpointing control.
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  S->RequireField(source_key_)
      ->SetMesh(S->GetMesh(source_domain_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(source_key_)->EnsureCompatibility(S);
      
}


} // namespace
} // namespace

