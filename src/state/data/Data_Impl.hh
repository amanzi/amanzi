/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_DATA_IMPL_HH_
#define AMANZI_DATA_IMPL_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "errors.hh"

#include "Data_Initializers.hh"
#include "Visualization.hh"
#include "Checkpoint.hh"

namespace Amanzi {

// underlying interface with type erasure
class Data_Intf {
 public:
  virtual ~Data_Intf() {}

  //  virtual std::unique_ptr<Data_Intf> Clone() const = 0;

  template <typename T>
  const T& Get() const;

  template <typename T>
  T& GetW();

  template <typename T>
  Teuchos::RCP<const T> GetPtr() const;

  template <typename T>
  Teuchos::RCP<T> GetPtrW();

  template <typename T>
  void SetPtr(Teuchos::RCP<T> t);

  template <typename T>
  void Set(const T& t);

  // virtual interface for ad-hoc polymorphism
  virtual void WriteVis(const Visualization& vis,
                        const Teuchos::ParameterList& attrs) const = 0;
  virtual void WriteCheckpoint(const Checkpoint& chkp,
                               const Teuchos::ParameterList& attrs) const = 0;
  virtual void ReadCheckpoint(const Checkpoint& chkp,
                              const Teuchos::ParameterList& attrs) = 0;
  virtual bool Initialize(Teuchos::ParameterList& plist,
                          const Teuchos::ParameterList& attrs) = 0;
};

// underlying implementation of type erasure
template <typename T>
class Data_Impl : public Data_Intf {
 public:
  Data_Impl() {}

  Data_Impl(Teuchos::RCP<T> t) : t_(std::move(t)) {}

  // Data_Impl(const Data_Impl<T>& other) :
  //     t_(Teuchos::rcp(new T(*other.t_))) {}

  Data_Impl(Data_Impl<T>&& other) noexcept : t_(std::move(other.t_)) {}

  template <typename... Ts>
  Data_Impl(Ts&&... ts) : t_(Teuchos::rcp(new T(std::forward<Ts...>(ts...))))
  {}

  // std::unique_ptr<Data_Intf> Clone() const override {
  //   return std::unique_ptr<Data_Impl<T> >(new Data_Impl<T>(*this));
  // }

  Data_Impl<T>& operator=(Data_Impl<T> other) = delete; // {

  const T& Get() const { return *t_; }
  T& GetW() { return *t_; }

  Teuchos::RCP<const T> GetPtr() const { return t_; }
  Teuchos::RCP<T> GetPtrW() { return t_; }

  void SetPtr(Teuchos::RCP<T> t) { t_ = std::move(t); }
  void Set(const T& t) { *t_ = t; }

  // virtual interface for ad-hoc polymorphism
  void WriteVis(const Visualization& vis,
                const Teuchos::ParameterList& attrs) const override
  {
    vis.Write(attrs, *t_);
  }

  void WriteCheckpoint(const Checkpoint& chkp,
                       const Teuchos::ParameterList& attrs) const override
  {
    chkp.Write(attrs, *t_);
  }

  void ReadCheckpoint(const Checkpoint& chkp,
                      const Teuchos::ParameterList& attrs) override
  {
    chkp.Read(attrs, *t_);
  }

  bool Initialize(Teuchos::ParameterList& plist,
                  const Teuchos::ParameterList& attrs) override
  {
    return Data_Initializers::Initialize(plist, attrs, *t_);
  }

 private:
  Teuchos::RCP<T> t_;
};

//
// Implementations of the underlying functions of the interface class, these
// cast to the correct type and provide access.
//

template <typename T>
const T&
Data_Intf::Get() const
{
  auto p = dynamic_cast<const Data_Impl<T>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type: \"" << typeid(T).name() << "\"";
    throw(msg);
  }
  return p->Get();
}

template <typename T>
T&
Data_Intf::GetW()
{
  auto p = dynamic_cast<Data_Impl<T>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type: \"" << typeid(T).name() << "\"";
    throw(msg);
  }
  return p->GetW();
}

template <typename T>
Teuchos::RCP<const T>
Data_Intf::GetPtr() const
{
  auto p = dynamic_cast<const Data_Impl<T>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type: \"" << typeid(T).name() << "\"";
    throw(msg);
  }
  return p->GetPtr();
}

template <typename T>
Teuchos::RCP<T>
Data_Intf::GetPtrW()
{
  auto p = dynamic_cast<Data_Impl<T>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type: \"" << typeid(T).name() << "\"";
    throw(msg);
  }
  return p->GetPtrW();
}

template <typename T>
void
Data_Intf::SetPtr(Teuchos::RCP<T> t)
{
  auto p = dynamic_cast<Data_Impl<T>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type: \"" << typeid(T).name() << "\"";
    throw(msg);
  }
  p->SetPtr(std::move(t));
}

template <typename T>
void
Data_Intf::Set(const T& t)
{
  auto p = dynamic_cast<Data_Impl<T>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << " data requested via incorrect type: \"" << typeid(T).name() << "\"";
    throw(msg);
  }
  p->Set(t);
}

} // namespace Amanzi

#endif
