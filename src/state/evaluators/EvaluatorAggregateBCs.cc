/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "Patch.hh"
#include "BCs.hh"
#include "DataStructuresHelpers.hh"

#include "EvaluatorAggregateBCs.hh"

namespace Amanzi {

const std::string EvaluatorAggregateBCs::name = "boundary condition aggregrator";

EvaluatorAggregateBCs::EvaluatorAggregateBCs(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondary(plist), inited_(false)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  auto deps = plist_->get<Teuchos::Array<std::string>>("dependencies");
  for (const auto& dep : deps) { dependencies_.insert(std::make_pair(dep, my_keys_[0].second)); }
}

void
EvaluatorAggregateBCs::EnsureCompatibility(State& S)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  auto my_keytag = my_keys_.front();
  auto& my_fac = S.Require<Operators::BCs, Operators::BCs_Factory>(
    my_keytag.first, my_keytag.second, my_keytag.first);

  for (auto& dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);

  // check plist for vis or checkpointing control
  EnsureCompatibility_Flags_(S);

  if (my_fac.mesh().get() && !inited_) {
    for (const auto& dep : dependencies_) {
      auto& eval = S.RequireEvaluator(dep.first, dep.second);
      auto& fac = S.Require<MultiPatch<double>, MultiPatchSpace>(dep.first, dep.second);
      fac.set_mesh(my_fac.mesh());
      fac.set_entity_kind(my_fac.entity_kind());
      eval.EnsureCompatibility(S);

      if (fac.flag_type == -1) {
        // also need a flag patch, but this is not in dependencies...
        auto& flag_fac = S.Require<MultiPatch<int>, MultiPatchSpace>(dep.first+"_flags", dep.second);
        flag_fac = fac;
      }
    }
    inited_ = true;
  }
}

void
EvaluatorAggregateBCs::Update_(State& S)
{
  auto& result = S.GetW<Operators::BCs>(my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);

  // overwrite with actual BCs
  {
    auto model = result.model();
    auto value = result.value();
    model->putScalar(0);
    value->putScalar(0.0);

    // set the default, 0 Neumann
    if (model->hasComponent("face")) {
      MultiVector_type_<int> model_bf(
        model->getMesh()->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false), 1);
      model_bf.putScalar(Operators::OPERATOR_BC_NEUMANN);
      model->getComponent("face", true)
        ->doExport(model_bf, model->getMesh()->getBoundaryFaceImporter(), Tpetra::INSERT);
    }

    // loop over dependencies and accumulate them
    for (const auto& dep : dependencies_) {
      const auto& i_bcs = S.Get<MultiPatch<double>>(dep.first, dep.second);
      copyMultiPatchToCompositeVector(i_bcs, to_string(result.kind()), *value, *model);

      if (i_bcs.space->flag_type == -1) {
        const auto& i_bcs_flags = S.Get<MultiPatch<int>>(dep.first+"_flags", dep.second);
        copyMultiPatchToCompositeVector(i_bcs_flags, to_string(result.kind()), *model);
      }
    }
  }

  // result.model()->print(std::cout);
  // result.value()->print(std::cout);
}

} // namespace Amanzi
