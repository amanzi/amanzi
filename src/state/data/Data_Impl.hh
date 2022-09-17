/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Basic object that may store heterogeneous data and put in State.
  The underlying object is stored as a Teuchos::RCP.

  Ownership: Note that this object TAKES OWNERSHIP of whatever it
  is created with. While initializing the data with a raw pointer
  is possible, this is a bad idea. Alternatively, you can pass in
  an RCP, which is then shared ownership.

  Implementation note -- Teuchos::RCP does NOT support move semantics,
  but we hope it will at some point.
*/

#ifndef AMANZI_STATE_DATA_IMPL_HH_
#define AMANZI_STATE_DATA_IMPL_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "errors.hh"

#include "Data_Helpers.hh"

namespace Amanzi {
namespace Impl {

// underlying interface with type erasure
class Data_Intf {
 public:
  virtual ~Data_Intf() {};

  // note, this clone is expected to return an empty Data object, so as not to
  // require a Clone method or copy constructor for the data itself -- that is
  // the job of the factory.
  virtual std::unique_ptr<Data_Intf> CloneEmpty() const = 0;

  template <typename T> const T& Get() const;
  template <typename T> T& GetW();

  template <typename T> Teuchos::RCP<const T> GetPtr() const;
  template <typename T> Teuchos::RCP<T> GetPtrW();

  template <typename T> void SetPtr(Teuchos::RCP<T> t);
  template <typename T> void Assign(const T& t);
  template <typename T> bool ValidType() const;

  // virtual interface for ad-hoc polymorphism
  virtual
  void WriteVis(const Visualization& vis, const Key& fieldname,
                const std::vector<std::string>* subfieldnames) const = 0;
  virtual
  void WriteCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                       const std::vector<std::string>* subfieldnames) const = 0;

  virtual
  bool ReadCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                      const std::vector<std::string>* subfieldnames) const = 0;

  virtual
  bool Initialize(Teuchos::ParameterList& plist, const Key& fieldname,
                  const std::vector<std::string>* subfieldnames) = 0;

  virtual
  void Assign(const Data_Intf& other) = 0;

  virtual
  const Data_Intf& operator=(const Data_Intf& other) = 0;
};


// underlying implementation of type erasure
template <typename T> class Data_Impl : public Data_Intf {
 public:
  Data_Impl() {};
  Data_Impl(Teuchos::RCP<T> t) : t_(std::move(t)) {};

  // Data_Impl(const Data_Impl<T>& other) :
  //     t_(Teuchos::rcp(new T(*other.t_))) {}

  Data_Impl(Data_Impl<T>&& other) noexcept : t_(std::move(other.t_)) {}

  template <typename... Ts>
  Data_Impl(Ts&&... ts)
      : t_(Teuchos::rcp(new T(std::forward<Ts...>(ts...)))) {}

  virtual std::unique_ptr<Data_Intf> CloneEmpty() const override {
    return std::unique_ptr<Data_Impl<T>>(new Data_Impl<T>());
  }

  virtual const Data_Intf& operator=(const Data_Intf& other) override {
    auto other_as_T = dynamic_cast<const Data_Impl<T>*>(&other);
    AMANZI_ASSERT(other_as_T != nullptr);
    t_ = other_as_T->t_;
    return *this;
  }

  const T& Get() const { return *t_; }
  T& GetW() { return *t_; }

  Teuchos::RCP<const T> GetPtr() const { return t_; }
  Teuchos::RCP<T> GetPtrW() { return t_; }

  void SetPtr(Teuchos::RCP<T> t) { t_ = std::move(t); }
  void Assign(const T& t) { *t_ = t; }

  // virtual interface for ad-hoc polymorphism
  virtual
  void WriteVis(const Visualization& vis, const Key& fieldname,
                const std::vector<std::string>* subfieldnames) const override {
    ::Amanzi::Helpers::WriteVis(vis, fieldname, subfieldnames, *t_);
  }

  virtual
  void WriteCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                       const std::vector<std::string>* subfieldnames) const override {
    ::Amanzi::Helpers::WriteCheckpoint(chkp, fieldname, subfieldnames, *t_);
  }

  virtual
  bool ReadCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                      const std::vector<std::string>* subfieldnames) const override {

    return ::Amanzi::Helpers::ReadCheckpoint(chkp, fieldname, subfieldnames, *t_);
  }

  virtual
  bool Initialize(Teuchos::ParameterList& plist, const Key& fieldname,
                  const std::vector<std::string>* subfieldnames) override {
    return ::Amanzi::Helpers::Initialize(plist, *t_, fieldname, subfieldnames);
  }

  virtual
  void Assign(const Data_Intf& other) override {
    return ::Amanzi::Helpers::Assign(*t_, other.Get<T>());
  }

 private:
  Teuchos::RCP<T> t_;
};


// Implementations of the underlying functions of the interface class, these
// cast to the correct type and provide access.
template <typename T> const T& Data_Intf::Get() const {
  auto p = dynamic_cast<const Data_Impl<T>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type (1): \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->Get();
}

template <typename T> T& Data_Intf::GetW() {
  auto p = dynamic_cast<Data_Impl<T>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type (2): \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->GetW();
}

template <typename T> Teuchos::RCP<const T> Data_Intf::GetPtr() const {
  auto p = dynamic_cast<const Data_Impl<T>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type (3): \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->GetPtr();
}

template <typename T> Teuchos::RCP<T> Data_Intf::GetPtrW() {
  auto p = dynamic_cast<Data_Impl<T>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type (4): \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->GetPtrW();
}

template <typename T> void Data_Intf::SetPtr(Teuchos::RCP<T> t) {
  auto p = dynamic_cast<Data_Impl<T>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type (5): \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  p->SetPtr(std::move(t));
}

template <typename T> void Data_Intf::Assign(const T& t) {
  auto p = dynamic_cast<Data_Impl<T>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type (6): \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  p->Assign(t);
}

template <typename T> bool Data_Intf::ValidType() const {
  auto p = dynamic_cast<const Data_Impl<T>* >(this);
  if (!p) return false;
  else return true;
}

} // namespace Impl
} // namespace Amanzi

#endif
