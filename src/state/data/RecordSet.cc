/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "RecordSet.hh"
#include "UniqueHelpers.hh"
#include "errors.hh"

namespace Amanzi {

// gather attributes for use with other functionality
Teuchos::ParameterList
RecordSet::attributes() const
{
  Teuchos::ParameterList attrs(vis_fieldname());
  attrs.set("units", units());
  attrs.set("location", location());
  attrs.set("subfieldnames", subfieldnames());
  return attrs;
}

// pass-throughs for other functionality
void
RecordSet::WriteVis(const Visualization& vis) const
{
  auto attrs = attributes();
  if (HasRecord("")) {
    records_.at("")->WriteVis(vis, attrs);
  } else if (HasRecord("next")) {
    records_.at("next")->WriteVis(vis, attrs);
  } else if (records_.size() > 0) {
    records_.begin()->second->WriteVis(vis, attrs);
  }
}
void
RecordSet::WriteCheckpoint(const Checkpoint& chkp) const
{
  if (HasRecord("")) {
    auto attrs = attributes();
    records_.at("")->WriteCheckpoint(chkp, attrs);
  }
}
void
RecordSet::ReadCheckpoint(const Checkpoint& chkp)
{
  if (HasRecord("")) {
    auto attrs = attributes();
    records_.at("")->ReadCheckpoint(chkp, attrs);
  }
}
bool
RecordSet::Initialize(Teuchos::ParameterList& plist)
{
  if (plist.isParameter("units")) set_units(plist.get<std::string>("units"));
  if (plist.isParameter("location"))
    set_location(AmanziMesh::entity_kind(plist.get<std::string>("location")));
  if (plist.isParameter("subfieldnames"))
    set_subfieldnames(plist.get<Teuchos::Array<std::string>>("subfieldnames"));

  bool init = false;
  if (HasRecord("")) {
    auto attrs = attributes();
    init = records_.at("")->Initialize(plist, attrs);
  }
  return init;
}

// Copy management
bool
RecordSet::HasRecord(const Key& tag) const
{
  return records_.count(tag) > 0;
}

Record&
RecordSet::GetRecord(const Key& tag)
{
  try {
    return *records_.at(tag);
  } catch (const std::out_of_range& e) {
    Errors::Message msg;
    msg << "Record: \"" << fieldname_ << "\" << does not have tag \"" << tag
        << "\"";
    throw(msg);
  }
}

const Record&
RecordSet::GetRecord(const Key& tag) const
{
  try {
    return *records_.at(tag);
  } catch (const std::out_of_range& e) {
    Errors::Message msg;
    msg << "Record: \"" << fieldname_ << "\" << does not have tag \"" << tag
        << "\"";
    throw(msg);
  }
}

// //  No data movement on pointer switch
// void RecordSet::SwitchCopies(const Key& tag1, const Key& tag2) {
//   GetRecord(tag1).swap(GetRecord(tag2));
// }

void
RecordSet::RequireRecord(const Key& tag, const Key& owner)
{
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

} // namespace Amanzi
