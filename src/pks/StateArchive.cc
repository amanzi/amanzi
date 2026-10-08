/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "EvaluatorPrimary.hh"

#include "StateArchive.hh"

namespace Amanzi {

/* ******************************************************************
* Deep copy of state fields.
****************************************************************** */
void
StateArchive::Add(const std::vector<std::string>& fields, const Tag& tag)
{
  tag_ = tag;
  for (const auto& name : fields) {
    fields_.emplace(name, S_->Get<CompositeVector>(name, tag));
  }
}


/* ******************************************************************
* Deep copy of state fields.
****************************************************************** */
void
StateArchive::Restore(const std::string& passwd)
{
  std::stringstream ss1, ss2;

  for (auto it = fields_.begin(); it != fields_.end(); ++it) {
    std::string owner(passwd);
    if (S_->HasEvaluator(it->first, tag_)) {
      auto type = S_->GetEvaluatorPtr(it->first, tag_)->get_type();
      if (type == EvaluatorType::SECONDARY || type == EvaluatorType::INDEPENDENT) owner = it->first;
    }
    S_->GetW<CompositeVector>(it->first, tag_, owner) = it->second;

    ss1 << "\"" << it->first << "\", ";

    if (S_->HasEvaluator(it->first, tag_)) {
      if (S_->GetEvaluatorPtr(it->first, tag_)->get_type() == EvaluatorType::PRIMARY) {
        Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(
          S_->GetEvaluatorPtr(it->first, tag_))
          ->SetChanged();

        ss2 << "\"" << it->first << "\", ";
      }
    }
  }

  if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "restored: " << ss1.str() << std::endl;
    *vo_->os() << "changed status of primaries: " << ss2.str() << std::endl;
  }
}


/* *******************************************************************
* Copy: Field (BASE) -> Field (prev_BASE)
* NOTE: boolean flag = false when the function is used in CommitStep()
******************************************************************* */
void
StateArchive::CopyFieldsToPrevFields(std::vector<std::string>& fields,
                                     const std::string& passwd,
                                     bool add)
{
  for (auto it = fields.begin(); it != fields.end(); ++it) {
    auto name = Keys::splitKey(*it);
    std::string prev = Keys::getKey(name.first, "prev_" + name.second);
    if (S_->HasRecord(prev, tag_)) {
      if (S_->HasEvaluator(*it, tag_)) {
        if (S_->GetEvaluatorPtr(*it, tag_)->get_type() == EvaluatorType::SECONDARY)
          S_->GetEvaluator(*it).Update(*S_, "archive");
      }
      if (add) fields_.emplace(prev, S_->Get<CompositeVector>(prev));
      S_->GetW<CompositeVector>(prev, tag_, passwd) = S_->Get<CompositeVector>(*it);

      // if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
      //   Teuchos::OSTab tab = vo_->getOSTab();
      //   *vo_->os() << "moved " << *it << " to " << prev << std::endl;
      // }
    }
  }
}


/* ******************************************************************
* Return a copy
****************************************************************** */
const CompositeVector&
StateArchive::get(const std::string& name)
{
  auto it = fields_.find(name);
  if (it != fields_.end() ) return it->second;
  AMANZI_ASSERT(false);
  // hide warnings
  return fields_.begin()->second;
}

} // namespace Amanzi
