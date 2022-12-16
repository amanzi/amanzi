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
  Record(){};
  Record(Key fieldname, Tag tag, Key owner = "");
  Record(const Record& other, const Tag* tag = nullptr);

  Record& operator=(const Record&) = default;
  Record& operator=(Record&&) = default;

  // access
  Key fieldname() const { return fieldname_; }
  bool initialized() const { return initialized_; }
  Key owner() const { return owner_; }
  Key vis_fieldname() const { return vis_key_; }
  bool io_checkpoint() const { return io_checkpoint_; }
  bool io_vis() const { return io_vis_; }

  // mutators
  void set_fieldname(const Key& fieldname) { fieldname_ = fieldname; }
  void set_initialized(bool initialized = true) { initialized_ = initialized; }
  void set_owner(const Key& owner) { owner_ = owner; }
  void set_vis_fieldname(const Key& key) { vis_key_ = key; }
  void set_io_checkpoint(bool io_checkpoint = true) { io_checkpoint_ = io_checkpoint; }
  void set_io_vis(bool io_vis = true) { io_vis_ = io_vis; }

  // pass-throughs for other functionality
  void
  WriteVis(const Visualization& vis, const std::vector<std::string>* subfieldnames = nullptr) const;
  void WriteCheckpoint(const Checkpoint& chkp,
                       const Tag& tag,
                       bool post_mortem = false,
                       const std::vector<std::string>* subfieldnames = nullptr) const;
  bool ReadCheckpoint(const Checkpoint& chkp,
                      const Tag& tag,
                      const std::vector<std::string>* subfieldnames = nullptr);
  bool Initialize(Teuchos::ParameterList& plist,
                  const std::vector<std::string>* subfieldnames = nullptr,
                  bool force = false);

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
  void Assign(const Key& owner, const T& t)
  {
    AssertOwnerOrDie(owner);
    try {
      data_.Assign<T>(t);
    } catch (const Errors::Message& msg) {
      Errors::Message new_msg;
      new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
      throw(new_msg);
    }
  }

  void Assign(const Record& other);
  void AssignPtr(const Record& other);

  // consistency checking
  void AssertOwnerOrDie(const Key& owner) const;

  template <typename T>
  bool ValidType() const
  {
    return data_.ValidType<T>();
  }

 private:
  Key fieldname_;
  Tag tag_;
  Key owner_;
  Key vis_key_;

  bool initialized_;
  bool io_checkpoint_;
  bool io_vis_;

  Impl::Data data_;

  friend class RecordSet;
};

} // namespace Amanzi

#endif
