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
StateArchive::Add(std::vector<std::string>& fields, const Tag& tag)
{
  tag_ = tag;
  for (const auto& name : fields) { fields_.emplace(name, S_->Get<CompositeVector>(name, tag)); }
}


/* ******************************************************************
* Deep copy of state fields.
****************************************************************** */
void
StateArchive::Restore(const std::string& passwd)
{
  for (auto it = fields_.begin(); it != fields_.end(); ++it) {
    if (S_->HasEvaluator(it->first, tag_)) {
      if (S_->GetEvaluatorPtr(it->first, tag_)->get_type() == EvaluatorType::SECONDARY)
        S_->GetW<CompositeVector>(it->first, tag_, it->first) = it->second;
    } else {
      S_->GetW<CompositeVector>(it->first, tag_, passwd) = it->second;
    }

    if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "reverted field \"" << it->first << "\"" << std::endl;
    }

    if (S_->HasEvaluator(it->first, tag_)) {
      if (S_->GetEvaluatorPtr(it->first, tag_)->get_type() == EvaluatorType::PRIMARY) {
        Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(
          S_->GetEvaluatorPtr(it->first, tag_))
          ->SetChanged();

        if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "changed status of primary field \"" << it->first << "\"" << std::endl;
        }
      }
    }
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
      if (add) fields_.emplace(prev, S_->Get<CompositeVector>(prev));
      S_->GetW<CompositeVector>(prev, tag_, passwd) = S_->Get<CompositeVector>(*it);
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
  if (it != fields_.end()) return it->second;

  AMANZI_ASSERT(false);
}

} // namespace Amanzi
