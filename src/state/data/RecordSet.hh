/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_STATE_DATA_RECORD_SET_HH_
#define AMANZI_STATE_DATA_RECORD_SET_HH_

#include <unordered_map>

#include "Teuchos_RCP.hpp"

#include "DataFactory.hh"
#include "Record.hh"
#include "StateDefs.hh"

namespace Amanzi {

class Visualization;


class RecordSet {
 private:
  using RecordMap = std::unordered_map<Key, std::unique_ptr<Record>>;

 public:
  // constructors
  RecordSet() {}
  RecordSet(Key fieldname)
    : fieldname_(fieldname), vis_fieldname_(fieldname), initialized_(false)
  {}

  // delete things that suggest this could be duplicated, etc
  RecordSet(const RecordSet& other) = delete;
  RecordSet(RecordSet&& other) = delete;
  RecordSet& operator=(const RecordSet&) = delete;
  RecordSet& operator=(RecordSet&&) = delete;

  // access
  const Key& fieldname() const { return fieldname_; }
  const Key& vis_fieldname() const { return vis_fieldname_; }
  bool initialized() const { return initialized_; }
  std::string units() const { return units_; }
  AmanziMesh::Entity_kind location() const { return location_; }
  const Teuchos::Array<std::string>& subfieldnames() const
  {
    return subfieldnames_;
  }
  Teuchos::ParameterList attributes() const;

  // mutate
  void set_fieldname(Key fieldname) { fieldname_ = std::move(fieldname); }
  void set_vis_fieldname(Key vis_fieldname)
  {
    vis_fieldname_ = std::move(vis_fieldname);
  }
  void set_initialized(bool initialized = true) { initialized_ = initialized; }
  void set_units(const std::string& units) { units_ = units; }
  void set_location(AmanziMesh::Entity_kind location) { location_ = location; }
  void set_subfieldnames(Teuchos::Array<std::string> subfieldnames)
  {
    subfieldnames_ = std::move(subfieldnames);
  }
  void set_io_checkpoint(bool io_checkpoint = true) {
    for (auto& e : records_) e.second->set_io_checkpoint(io_checkpoint);
  }
  void set_io_vis(bool io_vis = true) {
    for (auto& e : records_) e.second->set_io_vis(io_vis);
  }


  // pass-throughs for other functionality
  void WriteVis(const Visualization& vis, const Key& tag) const;
  void WriteCheckpoint(const Checkpoint& chkp) const;
  void ReadCheckpoint(const Checkpoint& chkp);
  bool Initialize(Teuchos::ParameterList& plist);

  // copy management
  const Record& GetRecord(const Key& tag) const;
  Record& GetRecord(const Key& tag);

  void RequireRecord(const Key& tag, const Key& owner);

  //  void SwitchCopies(const Key& tag1, const Key& tag2); // needs owner
  //  information?
  bool HasRecord(const Key& tag) const;
  bool DeleteRecord(const Key& tag);

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
      e.second->data_ = std::forward<Data>(factory_.Create());
    }
  }

  // Data setters/getters
  template <typename T>
  const T& Get(const Key& tag = "") const
  {
    return records_.at(tag)->Get<T>();
  }

  template <typename T>
  T& GetW(const Key& tag, const Key& owner)
  {
    return records_.at(tag)->GetW<T>(owner);
  }
  template <typename T>
  T& GetW(const Key& owner)
  {
    return GetW<T>("", owner);
  }

  template <typename T>
  Teuchos::RCP<const T> GetPtr(const Key& tag = "") const
  {
    return records_.at(tag)->GetPtr<T>();
  }

  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& tag, const Key& owner)
  {
    return records_.at(tag)->GetPtrW<T>(owner);
  }
  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& owner)
  {
    return GetPtrW<T>("", owner);
  }

  template <typename T>
  void SetPtr(const Key& tag, const Key& owner, const Teuchos::RCP<T>& t)
  {
    records_.at(tag)->SetPtr<T>(owner, t);
  }
  template <typename T>
  void SetPtr(const Key& owner, const Teuchos::RCP<T>& t)
  {
    SetPtr<T>("", owner, t);
  }

  template <typename T>
  void Set(const Key& tag, const Key& owner, const T& t)
  {
    records_.at(tag)->Set<T>(owner, t);
  }
  template <typename T>
  void Set(const Key& owner, const T& t)
  {
    Set<T>("", owner, t);
  }

  template <typename T, typename F>
  F& GetFactory()
  {
    return factory_.GetW<T, F>();
  }

  bool HasType() { return factory_.HasType(); }

  template <typename T, typename F>
  F& SetType()
  {
    if (!factory_.HasType()) { factory_ = dataFactory<T, F>(); }
    return GetFactory<T, F>();
  }

  template <typename T>
  void SetType()
  {
    if (!factory_.HasType()) { factory_ = dataFactory<T, NullFactory>(); }
    GetFactory<T, NullFactory>(); // checks valid type
  }

 private:
  Key fieldname_;
  Key vis_fieldname_;
  bool initialized_;
  std::string units_;
  AmanziMesh::Entity_kind location_;
  Teuchos::Array<std::string> subfieldnames_;

  DataFactory factory_;
  RecordMap records_;
  RecordMap derivatives_;
};

} // namespace Amanzi

#endif
