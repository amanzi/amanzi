/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Interface for a RecordSet. A field contains meta-data about 
  a data structure such as ownership, initialization, etc.
*/

#include "RecordSet.hh"
#include "UniqueHelpers.hh"
#include "errors.hh"

namespace Amanzi {

// pass-throughs for other functionality
void RecordSet::WriteVis(const Visualization& vis) const {
  for (auto& e : records_) {
    e.second->WriteVis(vis);
  }
}
void RecordSet::WriteCheckpoint(const Checkpoint& chkp) const {
  for (auto& e : records_) {
    e.second->WriteCheckpoint(chkp);
  }
}
void RecordSet::ReadCheckpoint(const Checkpoint& chkp) {
  for (auto& e : records_) {
    e.second->ReadCheckpoint(chkp);
  }
}
bool RecordSet::Initialize(Teuchos::ParameterList& plist) {
  bool init = false;
  for (auto& e : records_) {
    init |= e.second->Initialize(plist);
  }
  return init;
}

// Copy management
bool RecordSet::HasRecord(const Key& tag) const {
  return records_.count(tag) > 0;
}

Record& RecordSet::GetRecord(const Key& tag) {
  try {
    return *records_.at(tag);
  } catch (const std::out_of_range& e) {
    Errors::Message msg;
    msg << "Record: \"" << fieldname_ << "\" << does not have tag \"" << tag << "\"";
    throw(msg);
  }
}

const Record& RecordSet::GetRecord(const Key& tag) const {
  try {
    return *records_.at(tag);
  } catch (const std::out_of_range& e) {
    Errors::Message msg;
    msg << "Record: \"" << fieldname_ << "\" << does not have tag \"" << tag << "\"";
    throw(msg);
  }
}

void RecordSet::RequireRecord(const Key& tag, const Key& owner) {
  if (!HasRecord(tag)) {
    records_.emplace(tag, std::make_unique<Record>(fieldname(), owner));
    auto& r = records_.at(tag);
    if (!tag.empty()) {
      r->set_vis_fieldname(vis_fieldname() + std::string("_") + tag);
    }
  } else if (!owner.empty()) {
    auto& r = records_.at(tag);
    if (r->owner().empty()) {
      r->set_owner(owner);
    } else {
      r->AssertOwnerOrDie(owner);
    }
  }
}

}  // namespace Amanzi
