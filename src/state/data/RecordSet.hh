/* -------------------------------------------------------------------------
  Arcos

  License: see $ATS_DIR/COPYRIGHT
  Author: Ethan Coon

  Interface for a Record.  A record contains meta-data about a data
  structure and that data structure, encapsulating things like vis,
  checkpointing, ownership, initialization, etc.
  ------------------------------------------------------------------------- */

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
      : fieldname_(fieldname), vis_fieldname_(fieldname), initialized_(false) {}

  // delete things that suggest this could be duplicated, etc
  RecordSet(const RecordSet &other) = delete;
  RecordSet(RecordSet &&other) = delete;
  RecordSet &operator=(const RecordSet &) = delete;
  RecordSet &operator=(RecordSet &&) = delete;

  // access
  const Key &fieldname() const { return fieldname_; }
  const Key &vis_fieldname() const { return vis_fieldname_; }
  bool initialized() const { return initialized_; }
  Units units() const { return units_; }

  // mutate
  void set_fieldname(Key fieldname) { fieldname_ = std::move(fieldname); }
  void set_vis_fieldname(Key vis_fieldname) {
    vis_fieldname_ = std::move(vis_fieldname);
  }
  void set_initialized(bool initialized = true) { initialized_ = initialized; }
  void set_units(Units units) { units = std::move(units); }

  // pass-throughs for other functionality
  void WriteVis(const Visualization &vis) const;
  void WriteCheckpoint(const Checkpoint &chkp) const;
  void ReadCheckpoint(const Checkpoint &chkp);
  bool Initialize(Teuchos::ParameterList &plist);

  // copy management
  const Record &GetRecord(const Key &tag) const;
  Record &GetRecord(const Key &tag);

  void RequireRecord(const Key &tag, const Key &owner);
  void RequireDerivativeRecord(const Key &tag, const Key &wrt_key,
                               const Key &wrt_tag, const Key &owner);

  //  void SwitchCopies(const Key& tag1, const Key& tag2); // needs owner
  //  information?
  bool HasRecord(const Key &tag) const;
  bool HasDerivativeRecord(const Key &tag, const Key &wrt_key,
                           const Key &wrt_tag) const;
  bool DeleteRecord(const Key &tag);

  // Iterate over tags
  typedef RecordMap::const_iterator tag_iterator;
  tag_iterator begin() const { return records_.begin(); }
  tag_iterator end() const { return records_.end(); }
  RecordMap::size_type count() { return records_.size(); }
  RecordMap::size_type size() { return records_.size(); }

  // Data creation
  void CreateData() {
    for (auto &e : records_) {
      e.second->data_ = std::forward<Data>(factory_.Create());
    }
  }

  // Data setters/getters
  template <typename T> const T &Get(const Key &tag = "") const {
    return records_.at(tag)->Get<T>();
  }

  template <typename T> T &GetW(const Key &tag, const Key &owner) {
    return records_.at(tag)->GetW<T>(owner);
  }
  template <typename T> T &GetW(const Key &owner) { return GetW<T>("", owner); }

  template <typename T>
  Teuchos::RCP<const T> GetPtr(const Key &tag = "") const {
    return records_.at(tag)->GetPtr<T>();
  }

  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key &tag, const Key &owner) {
    return records_.at(tag)->GetPtrW<T>(owner);
  }
  template <typename T> Teuchos::RCP<T> GetPtrW(const Key &owner) {
    return GetPtrW<T>("", owner);
  }

  template <typename T>
  void SetPtr(const Key &tag, const Key &owner, const Teuchos::RCP<T> &t) {
    records_.at(tag)->SetPtr<T>(owner, t);
  }
  template <typename T>
  void SetPtr(const Key &owner, const Teuchos::RCP<T> &t) {
    SetPtr<T>("", owner, t);
  }

  template <typename T> void Set(const Key &tag, const Key &owner, const T &t) {
    records_.at(tag)->Set<T>(owner, t);
  }
  template <typename T> void Set(const Key &owner, const T &t) {
    Set<T>("", owner, t);
  }

  template <typename T>
  const T &GetDerivative(const Key &wrt_key, const Key &wrt_tag) const {
    Key deriv = std::string{"d_d"} + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->Get<T>();
  }

  template <typename T>
  const T &GetDerivative(const Key &tag, const Key &wrt_key,
                         const Key &wrt_tag) const {
    Key deriv = std::string{"d"} + tag + "_d" + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->Get<T>();
  }

  template <typename T>
  T &GetDerivativeW(const Key &wrt_key, const Key &wrt_tag, const Key &owner) {
    Key deriv = std::string{"d_d"} + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->GetW<T>(owner);
  }

  template <typename T>
  T &GetDerivativeW(const Key &tag, const Key &wrt_key, const Key &wrt_tag,
                    const Key &owner) {
    Key deriv = std::string{"d"} + tag + "_d" + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->GetW<T>(owner);
  }

  template <typename T>
  Teuchos::RCP<const T> GetDerivativePtr(const Key &wrt_key,
                                         const Key &wrt_tag) const {
    Key deriv = std::string{"d_d"} + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->GetPtr<T>();
  }

  template <typename T>
  Teuchos::RCP<const T> GetDerivativePtr(const Key &tag, const Key &wrt_key,
                                         const Key &wrt_tag) const {
    Key deriv = std::string{"d"} + tag + "_d" + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->GetPtr<T>();
  }

  template <typename T>
  Teuchos::RCP<T> GetDerivativePtrW(const Key &wrt_key, const Key &wrt_tag,
                                    const Key &owner) {
    Key deriv = std::string{"d_d"} + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->GetPtrW<T>(owner);
  }

  template <typename T>
  Teuchos::RCP<T> GetDerivativePtrW(const Key &tag, const Key &wrt_key,
                                    const Key &wrt_tag, const Key &owner) {
    Key deriv = std::string{"d"} + tag + "_d" + wrt_key + ":" + wrt_tag;
    return derivatives_.at(deriv)->GetPtrW<T>(owner);
  }

  template <typename T, typename F> F &GetFactory() {
    return factory_.GetW<T, F>();
  }

  bool HasType() { return factory_.HasType(); }

  template <typename T, typename F> F &SetType() {
    if (!factory_.HasType()) {
      factory_ = dataFactory<T, F>();
    }
    return GetFactory<T, F>();
  }

  template <typename T> void SetType() {
    if (!factory_.HasType()) {
      factory_ = dataFactory<T, NullFactory>();
    }
    GetFactory<T, NullFactory>(); // checks valid type
  }

private:
  Key fieldname_;
  Key vis_fieldname_;
  bool initialized_;
  Units units_;

  DataFactory factory_;
  RecordMap records_;
  RecordMap derivatives_;
};

} // namespace Amanzi

#endif
