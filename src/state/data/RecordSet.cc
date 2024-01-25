/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Interface for a RecordSet. A field contains meta-data about
  a data structure such as ownership, initialization, etc.
*/

#include <memory>

#include "RecordSet.hh"
#include "errors.hh"

namespace Amanzi {

// gather attributes for use with other functionality
Teuchos::ParameterList
RecordSet::attributes(Tag const* const tag = nullptr) const
{
  std::string name;
  if (tag) name = Keys::getKey(vis_fieldname(), *tag);
  else name = vis_fieldname();

  Teuchos::ParameterList attrs(name);
  // attrs.set("units", units());
  attrs.set("location", location());
  auto names = subfieldnames();
  if (names) attrs.set("subfieldnames", names);
  return attrs;
}


// pass-throughs for other functionality
void
RecordSet::WriteVis(const Visualization& vis, Tag const* const tag) const
{
  if (tag) {
    if (HasRecord(*tag)) {
      auto attrs = attributes();
      GetRecord(*tag).WriteVis(vis, attrs);
    }
  } else {
    for (auto& e : records_) {
      auto attrs = attributes(&e.first);
      e.second->WriteVis(vis, attrs);
    }
  }
}

void
RecordSet::WriteCheckpoint(const Checkpoint& chkp, bool post_mortem) const
{
  for (auto& e : records_) {
    auto attrs = attributes(&e.first);
    e.second->WriteCheckpoint(chkp, attrs, post_mortem);
  }
}

void
RecordSet::ReadCheckpoint(const Checkpoint& chkp, Tag const* const tag)
{
  if (tag && HasRecord(*tag)) {
    auto attrs = attributes(tag);
    GetRecord(*tag).ReadCheckpoint(chkp, attrs);
  } else {
    for (auto& e : records_) {
      auto attrs = attributes(&e.first);
      e.second->ReadCheckpoint(chkp, attrs);
    }
  }
}

bool
RecordSet::Initialize(Teuchos::ParameterList& plist, bool force)
{
  bool init = false;

  if (plist.isParameter("units")) set_units(plist.get<std::string>("units"));
  if (plist.isParameter("location"))
    set_location(AmanziMesh::createEntityKind(plist.get<std::string>("location")));
  if (plist.isParameter("subfieldnames")) {
    set_subfieldnames(plist.get<Teuchos::Array<std::string>>("subfieldnames"));
  } else {
    auto names = subfieldnames();
    if (names)
      plist.set<Teuchos::Array<std::string>>("subfieldnames", *names);
  }

  for (auto& e : records_) { init |= e.second->Initialize(plist, force); }
  return init;
}

void
RecordSet::Assign(const Tag& dest, const Tag& source)
{
  GetRecord(dest).Assign(GetRecord(source));
}

void
RecordSet::AssignPtr(const Tag& alias, const Tag& target)
{
  GetRecord(alias).AssignPtr(GetRecord(target));
}


// Copy management
bool
RecordSet::HasRecord(const Tag& tag) const
{
  return records_.count(tag) > 0;
}

Record&
RecordSet::GetRecord(const Tag& tag)
{
  try {
    return *records_.at(tag);
  } catch (const std::out_of_range& e) {
    Errors::Message msg;
    msg << "Record: \"" << fieldname_ << "\" << does not have tag \"" << tag.get() << "\"";
    throw(msg);
  }
}

const Record&
RecordSet::GetRecord(const Tag& tag) const
{
  try {
    return *records_.at(tag);
  } catch (const std::out_of_range& e) {
    Errors::Message msg;
    msg << "Record: \"" << fieldname_ << "\" << does not have tag \"" << tag.get() << "\"";
    throw(msg);
  }
}


void
RecordSet::AliasRecord(const Tag& target, const Tag& alias)
{
  records_[alias] = std::make_shared<Record>(*records_[target], &alias);
  aliases_[alias] = target;
}


Record&
RecordSet::RequireRecord(const Tag& tag, const Key& owner, bool alias_ok)
{
  // Check if this should be an aliased record
  if (!HasRecord(tag) && alias_ok && tag == Tags::NEXT) {
    for (const auto& other_tag : records_) {
      if (other_tag.first != Tags::NEXT && Keys::in(other_tag.first.get(), "next")) {
        // alias!
        AliasRecord(other_tag.first, Tags::NEXT);
        return *other_tag.second;
      }
    }
  }

  // otherwise create the record
  if (!HasRecord(tag)) {
    records_.emplace(tag, std::make_shared<Record>(fieldname(), tag, owner));
    auto& r = records_.at(tag);
    // if (!tag.get().empty()) {
    //   r->set_vis_fieldname(Keys::getKey(vis_fieldname(), tag));
    // }
    return *r;
  } else {
    auto& r = records_.at(tag);
    if (!owner.empty()) {
      if (r->owner().empty()) {
        r->set_owner(owner);
      } else {
        r->AssertOwnerOrDie(owner);
      }
    }
    return *r;
  }
}


bool
RecordSet::isInitialized(Tag& failed)
{
  for (auto& r : records_) {
    if (!r.second->initialized()) {
      ;
      failed = r.first;
      return false;
    }
  }
  return true;
}

} // namespace Amanzi
