/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "BCs.hh"
#include "Evaluator_BCs.hh"
#include "Patch.hh"

namespace Amanzi {

Evaluator_BCs::Evaluator_BCs(Teuchos::ParameterList& plist)
    : EvaluatorSecondary(plist),
      inited_(false)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  for (auto sublist : plist.sublist("boundary functions")) {
    std::string sublist_name = Keys::cleanPListName(sublist.first);
    dependencies_.push_back(std::make_pair(sublist_name, my_keys_[0].second));
  }
}

void
Evaluator_BCs::EnsureCompatibility(State& S)
{
  auto& my_fac = S.Require<Operators::BCs, Operators::BCs_Factory>(
    my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);

  for (auto& dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);

  // check plist for vis or checkpointing control
  EnsureCompatibility_Flags_(S);

  if (my_fac.mesh().get() && !inited_) {
    for (const auto& dep : dependencies_) {
      auto& eval = S.RequireEvaluator(dep.first, dep.second);
      auto& fac = S.Require<MultiPatch, MultiPatchSpace>(dep.first, dep.second);
      fac.set_mesh(my_fac.mesh());
      fac.flag_entity = my_fac.entity_kind();
      eval.EnsureCompatibility(S);

    }
    inited_ = true;
  }
}

void
Evaluator_BCs::Update_(State& S)
{
  auto& result = S.GetW<Operators::BCs>(
    my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);

  // overwrite with actual BCs
  {
    auto model = result.model();
    auto value = result.value();
    model->putScalar(0);
    value->putScalar(0.0);

    // set the default, 0 Neumann
    if (model->HasComponent("face")) {
      MultiVector_type_<int> model_bf(model->Mesh()->exterior_face_map(true), 1);
      model_bf.putScalar(Operators::OPERATOR_BC_NEUMANN);
      model->GetComponent("face", true)->doExport(model_bf, *model->Mesh()->exterior_face_importer(), Tpetra::INSERT);
    }

    // loop over dependencies and accumulate them
    for (const auto& dep : dependencies_) {
      const auto& i_bcs = S.Get<MultiPatch>(dep.first, dep.second);
      multiPatchToCompositeVector(i_bcs, AmanziMesh::entity_kind_string(result.kind()), *value, *model);
    }
  }
}

} // namespace Amanzi
