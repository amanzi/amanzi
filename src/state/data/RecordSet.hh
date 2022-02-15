/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

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
  using RecordMap = std::unordered_map<Tag, std::unique_ptr<Record> >;

 public:
  // constructors
  RecordSet() {}
  RecordSet(Key fieldname)
      : fieldname_(fieldname),
        vis_fieldname_(fieldname) {}

  // delete things that suggest this could be duplicated, etc
  RecordSet(const RecordSet& other) = delete;
  RecordSet(RecordSet&& other) = delete;
  RecordSet& operator=(const RecordSet&) = delete;
  RecordSet& operator=(RecordSet&&) = delete;

  // access
  const Key& fieldname() const { return fieldname_; }
  const Key& vis_fieldname() const { return vis_fieldname_; }
  Utils::Units units() const { return units_; }

  // mutate
  void set_fieldname(Key fieldname) { fieldname_ = std::move(fieldname); }
  void set_vis_fieldname(Key vis_fieldname) {
    vis_fieldname_ = std::move(vis_fieldname);
  }
  void set_units(Utils::Units units) { units_ = std::move(units); }

  // pass-throughs for other functionality
  void WriteVis(const Visualization& vis) const;
  void WriteCheckpoint(const Checkpoint& chkp) const;
  void ReadCheckpoint(const Checkpoint& chkp);
  bool Initialize(Teuchos::ParameterList& plist);
  void Assign(const Tag& dest, const Tag& source);

  // copy management
  const Record& GetRecord(const Tag& tag) const;
  Record& GetRecord(const Tag& tag);

  Record& RequireRecord(const Tag& tag, const Key& owner);

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
  void CreateData() {
    for (auto& e : records_) {
      e.second->data_ = std::forward<Data>(factory_.Create());
    }
  }

  // Data setters/getters
  template <typename T> const T& Get(const Tag& tag = Tags::DEFAULT) const {
    return GetRecord(tag).Get<T>();
  }

  template <typename T> T& GetW(const Tag& tag, const Key& owner) {
    return GetRecord(tag).GetW<T>(owner);
  }
  template <typename T> T& GetW(const Key& owner) { return GetW<T>(Tags::DEFAULT, owner); }

  template <typename T>
  Teuchos::RCP<const T> GetPtr(const Tag& tag = Tags::DEFAULT) const {
    return GetRecord(tag).GetPtr<T>();
  }

  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Tag& tag, const Key& owner) {
    return GetRecord(tag).GetPtrW<T>(owner);
  }
  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& owner) {
    return GetPtrW<T>(Tags::DEFAULT, owner);
  }

  template <typename T>
  void SetPtr(const Tag& tag, const Key& owner, const Teuchos::RCP<T>& t) {
    GetRecord(tag).SetPtr<T>(owner, t);
  }
  template <typename T>
  void SetPtr(const Key& owner, const Teuchos::RCP<T>& t) {
    SetPtr<T>(Tags::DEFAULT, owner, t);
  }

  template <typename T>
  void Assign(const Tag& tag, const Key& owner, const T& t) {
    GetRecord(tag).Assign<T>(owner, t);
  }
  template <typename T>
  void Assign(const Key& owner, const T& t) {
    Assign<T>(Tags::DEFAULT, owner, t);
  }

  template <typename T, typename F>
  F& GetFactory() {
    return factory_.GetW<T, F>();
  }

  bool HasType() { return factory_.HasType(); }

  template <typename T, typename F>
  F& SetType() {
    if (!factory_.HasType()) {
      factory_ = dataFactory<T, F>();
    }
    return GetFactory<T, F>();
  }

  template <typename T>
  void SetType() {
    if (!factory_.HasType()) {
      factory_ = dataFactory<T, NullFactory>();
    }
    GetFactory<T, NullFactory>();  // checks valid type
  }

  template <typename T, typename F>
  bool ValidType() {
    return factory_.ValidType<T,F>();
  }

  template <typename T>
  bool ValidType() {
    return factory_.ValidType<T>();
  }

  // initialization of set is the collective (AND) operation
  bool isInitialized(Tag& failed);
  void initializeTags(bool initialized = true) {
    for (auto& r : records_) r.second->set_initialized(initialized);
  }

 private:
  Key fieldname_;
  Key vis_fieldname_;
  Utils::Units units_;

  DataFactory factory_;
  RecordMap records_;
};

} // namespace Amanzi

#endif
