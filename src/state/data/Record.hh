/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_STATE_DATA_RECORD_HH_
#define AMANZI_STATE_DATA_RECORD_HH_

#include <unordered_map>

#include "Teuchos_RCP.hpp"

#include "Data.hh"
#include "StateDefs.hh"

namespace Amanzi {

class Visualization;

class Record {
 public:
  Record() {}
  Record(Key fieldname, Key owner = "");

  Record(const Record& other);
  Record(Record&& other) = default;
  Record(const Record& other, Key fieldname) : Record(other)
  {
    set_fieldname(fieldname);
  }
  Record(const Record& other, Key fieldname, Key owner)
    : Record(other, fieldname)
  {
    set_owner(owner);
  }

  Record& operator=(const Record&) = default;
  Record& operator=(Record&&) = default;

  // access
  Key fieldname() const { return fieldname_; }
  bool initialized() const { return initialized_; }
  Key owner() const { return owner_; }
  Key vis_fieldname() const { return vis_key_; }
  bool io_checkpoint() const { return io_checkpoint_; }
  bool io_vis() const { return io_vis_; }
  void attributes(Teuchos::ParameterList& plist) const
  {
    plist.setName(vis_fieldname());
  }

  // mutators
  void set_fieldname(Key fieldname) { fieldname_ = std::move(fieldname); }
  void set_initialized(bool initialized = true) { initialized_ = initialized; }
  void set_owner(Key owner) { owner_ = std::move(owner); }
  void set_vis_fieldname(Key key) { vis_key_ = std::move(key); }
  void set_io_checkpoint(bool io_checkpoint = true)
  {
    io_checkpoint_ = std::move(io_checkpoint);
  }
  void set_io_vis(bool io_vis = true) { io_vis_ = std::move(io_vis); }

  // pass-throughs for other functionality
  void WriteVis(const Visualization& vis, Teuchos::ParameterList attrs) const;
  void
  WriteCheckpoint(const Checkpoint& chkp, Teuchos::ParameterList attrs) const;
  void ReadCheckpoint(const Checkpoint& chkp, Teuchos::ParameterList attrs);
  bool Initialize(Teuchos::ParameterList& plist, Teuchos::ParameterList attrs);

  // Data setters/getters
  template <typename T>
  const T& Get() const
  {
    try {
      return data_.Get<T>();
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  template <typename T>
  T& GetW(const Key& owner)
  {
    AssertOwnerOrDie(owner);
    try {
      return data_.GetW<T>();
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  template <typename T>
  Teuchos::RCP<const T> GetPtr() const
  {
    try {
      return data_.GetPtr<T>();
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& owner)
  {
    AssertOwnerOrDie(owner);
    try {
      return data_.GetPtrW<T>();
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  template <typename T>
  void SetPtr(const Key& owner, const Teuchos::RCP<T>& t)
  {
    AssertOwnerOrDie(owner);
    try {
      data_.SetPtr(t);
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  template <typename T>
  void Set(const Key& owner, const T& t)
  {
    AssertOwnerOrDie(owner);
    try {
      data_.Set(t);
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  // consistency checking
  void AssertOwnerOrDie(const Key& owner) const;

 private:
  Key fieldname_;
  Key owner_;
  Key vis_key_;
  std::vector<std::string> subfieldnames_;
  AmanziMesh::Entity_kind location_;

  bool initialized_;
  bool io_checkpoint_;
  bool io_vis_;

  Units units_;

  Data data_;

  friend class RecordSet;
};

} // namespace Amanzi

#endif
