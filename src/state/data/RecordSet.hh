/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Interface for a Record.  A record contains meta-data about a data
  structure and that data structure, encapsulating things like vis,
  checkpointing, ownership, initialization, etc.
*/

#ifndef AMANZI_STATE_DATA_RECORD_SET_HH_
#define AMANZI_STATE_DATA_RECORD_SET_HH_

#include <unordered_map>

#include "Teuchos_RCP.hpp"

#include "DataFactory.hh"
#include "Record.hh"
#include "StateDefs.hh"
#include "Tag.hh"

namespace Amanzi {

class Visualization;

class RecordSet {
 private:
  // this must be std::map because we loop over these for checkpoint/restart, which
  // requires a known (alphabetic) order across ranks.
  using RecordMap = std::map<Tag, std::shared_ptr<Record>>;

 public:
  // constructors
  RecordSet() {}
  RecordSet(Key fieldname) : fieldname_(fieldname), vis_fieldname_(fieldname) {}

  // delete things that suggest this could be duplicated, etc
  RecordSet(const RecordSet& other) = delete;
  RecordSet(RecordSet&& other) = delete;
  RecordSet& operator=(const RecordSet&) = delete;
  RecordSet& operator=(RecordSet&&) = delete;

  // access
  const Key& fieldname() const { return fieldname_; }
  const Key& vis_fieldname() const { return vis_fieldname_; }
  Utils::Units units() const { return units_; }
  std::vector<std::string> const* subfieldnames() const { return subfieldnames_.get(); }

  // mutate
  void set_fieldname(const Key& fieldname) { fieldname_ = fieldname; }
  void set_vis_fieldname(const Key& vis_fieldname) { vis_fieldname_ = vis_fieldname; }
  void set_units(const Utils::Units& units) { units_ = units; }
  void set_subfieldnames(const std::vector<std::string>& subfieldnames)
  {
    if (!subfieldnames_)
      subfieldnames_ = std::make_unique<std::vector<std::string>>(subfieldnames);
    else
      *subfieldnames_ = subfieldnames;
  }

  // pass-throughs for other functionality
  void WriteVis(const Visualization& vis, Tag const* const = nullptr) const;
  void WriteCheckpoint(const Checkpoint& chkp, bool post_mortem = false) const;
  void ReadCheckpoint(const Checkpoint& chkp);
  bool Initialize(Teuchos::ParameterList& plist, bool force = false);
  void Assign(const Tag& dest, const Tag& source);
  void AssignPtr(const Tag& dest, const Tag& source);

  // copy management
  const Record& GetRecord(const Tag& tag) const;
  Record& GetRecord(const Tag& tag);
  void AliasRecord(const Tag& target, const Tag& alias);

  Record& RequireRecord(const Tag& tag, const Key& owner, bool alias_ok = true);

  //  void SwitchCopies(const Tag& tag1, const Tag& tag2); // needs owner
  //  information?
  bool HasRecord(const Tag& tag) const;
  bool DeleteRecord(const Tag& tag);

  // Iterate over tags
  typedef RecordMap::const_iterator tag_iterator;
  tag_iterator begin() const { return records_.begin(); }
  tag_iterator end() const { return records_.end(); }
  RecordMap::size_type count() { return records_.size(); }
  RecordMap::size_type size() { return records_.size(); }

  // Data creation
  void CreateData()
  {
    for (auto& e : records_) {
      if (!aliases_.count(e.first)) {
        e.second->data_ = std::forward<Impl::Data>(factory_.Create());
      }
    }
    for (auto& pair : aliases_) { records_[pair.first]->AssignPtr(*records_[pair.second]); }
  }

  // Data setters/getters
  template <typename T>
  const T& Get(const Tag& tag = Tags::DEFAULT) const
  {
    return GetRecord(tag).Get<T>();
  }

  template <typename T>
  T& GetW(const Tag& tag, const Key& owner)
  {
    return GetRecord(tag).GetW<T>(owner);
  }
  template <typename T>
  T& GetW(const Key& owner)
  {
    return GetW<T>(Tags::DEFAULT, owner);
  }

  template <typename T>
  Teuchos::RCP<const T> GetPtr(const Tag& tag = Tags::DEFAULT) const
  {
    return GetRecord(tag).GetPtr<T>();
  }

  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Tag& tag, const Key& owner)
  {
    return GetRecord(tag).GetPtrW<T>(owner);
  }
  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& owner)
  {
    return GetPtrW<T>(Tags::DEFAULT, owner);
  }

  template <typename T>
  void SetPtr(const Tag& tag, const Key& owner, const Teuchos::RCP<T>& t)
  {
    GetRecord(tag).SetPtr<T>(owner, t);
  }
  template <typename T>
  void SetPtr(const Key& owner, const Teuchos::RCP<T>& t)
  {
    SetPtr<T>(Tags::DEFAULT, owner, t);
  }

  template <typename T>
  void Assign(const Tag& tag, const Key& owner, const T& t)
  {
    GetRecord(tag).Assign<T>(owner, t);
  }
  template <typename T>
  void Assign(const Key& owner, const T& t)
  {
    Assign<T>(Tags::DEFAULT, owner, t);
  }

  template <typename T, typename F>
  F& GetFactory()
  {
    return factory_.GetW<T, F>();
  }

  bool HasType() { return factory_.HasType(); }

  template <typename T, typename F>
  F& SetType(const F& f)
  {
    if (!factory_.HasType()) {
      factory_ = Impl::dataFactory<T, F>(f);
    } else {
      if (!Helpers::Equivalent(f, GetFactory<T, F>())) {
        Errors::Message msg;
        msg << "Factory required for field \"" << fieldname_
            << "\" differs from previous requirement call.";
        Exceptions::amanzi_throw(msg);
      }
    }
    return GetFactory<T, F>();
  }

  template <typename T, typename F>
  F& SetType()
  {
    if (!factory_.HasType()) { factory_ = Impl::dataFactory<T, F>(); }
    return GetFactory<T, F>();
  }

  template <typename T>
  void SetType()
  {
    if (!factory_.HasType()) { factory_ = Impl::dataFactory<T, NullFactory>(); }
    GetFactory<T, NullFactory>(); // checks valid type
  }

  template <typename T, typename F>
  bool ValidType() const
  {
    return factory_.ValidType<T, F>();
  }

  template <typename T>
  bool ValidType() const
  {
    return factory_.ValidType<T>();
  }

  // initialization of set is the collective (AND) operation
  bool isInitialized(Tag& failed);
  void initializeTags(bool initialized = true)
  {
    for (auto& r : records_) r.second->set_initialized(initialized);
  }

 private:
  Key fieldname_;
  Key vis_fieldname_;
  Utils::Units units_;
  std::unique_ptr<std::vector<std::string>> subfieldnames_;

  std::map<Tag, Tag> aliases_;
  Impl::DataFactory factory_;
  RecordMap records_;
};

} // namespace Amanzi

#endif
