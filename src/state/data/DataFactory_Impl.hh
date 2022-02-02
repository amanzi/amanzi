/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  Basic factory wrapper that may store heterogeneous data and
  create heterogeneous Data objects to put into State.

  This shares the same design pattern as Data to do type erasure.

  Implementation details of the underlying classes that do type erasure.
  If I were smarter, these could probably be made private members of
  the public class, but then the non-member constructors would need to
  be friends, etc.
*/

#ifndef AMANZI_STATE_STOREABLE_FACTORY_IMPL_HH_
#define AMANZI_STATE_STOREABLE_FACTORY_IMPL_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "StateDefs.hh"
#include "dbc.hh"

namespace Amanzi {

// thing factory base class
class DataFactory_Intf {
 public:
  virtual ~DataFactory_Intf() {};

  virtual std::unique_ptr<DataFactory_Intf> Clone() const = 0;

  template <typename T, typename F> bool ValidType() const;

  virtual void Create(Data& t) const = 0;
  virtual Data Create() const = 0;

  template <typename T, typename F> const F& Get() const;

  template <typename T, typename F> F& GetW();
};


// thing factory implementation
template <typename T, typename F>
class DataFactory_Impl : public DataFactory_Intf {
public:
  DataFactory_Impl() : f_(std::move(std::unique_ptr<F>(new F()))) {}

  DataFactory_Impl(const DataFactory_Impl& other)
      : f_(std::move(std::unique_ptr<F>(new F(*other.f_)))) {}

  std::unique_ptr<DataFactory_Intf> Clone() const override {
    return std::unique_ptr<DataFactory_Impl<T, F> >(
        new DataFactory_Impl<T, F>(*this));
  }

  void Create(Data& t) const override { t.SetPtr<T>(f_->Create()); }

  Data Create() const override { return data<T>(f_->Create()); }

  const F& Get() const { return *f_; }

  F& GetW() { return *f_; }

 private:
  std::unique_ptr<F> f_;
};


// partial specialization for NullFactory
template <typename T>
class DataFactory_Impl<T, NullFactory> : public DataFactory_Intf {
public:
  DataFactory_Impl()
      : f_(std::unique_ptr<NullFactory>(new NullFactory())) {};

  DataFactory_Impl(const DataFactory_Impl& other)
      : f_(std::unique_ptr<NullFactory>(new NullFactory(*other.f_))) {};

  std::unique_ptr<DataFactory_Intf> Clone() const override {
    return std::unique_ptr<DataFactory_Impl<T, NullFactory>>(
        new DataFactory_Impl<T, NullFactory>(*this));
  }

  // factory Create specialization on NullFactory factory (i.e. null factory)
  void Create(Data& t) const override { t.SetPtr<T>(Teuchos::rcp(new T())); }

  Data Create() const override { return data<T>(Teuchos::rcp(new T())); }

  const NullFactory& Get() const { return *f_; }

  NullFactory& GetW() { return *f_; }

 private:
  std::unique_ptr<NullFactory> f_;
};


// thing factory base class implementation
template <typename T, typename F> bool DataFactory_Intf::ValidType() const {
  auto p = dynamic_cast<const DataFactory_Impl<T, F>* >(this);
  if (!p) {
    return false;
  } else {
    return true;
  }
}

template <typename T, typename F> const F& DataFactory_Intf::Get() const {
  auto p = dynamic_cast<DataFactory_Impl<T, F>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << "factory requested via incorrect type: \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->Get();
}

template <typename T, typename F> F& DataFactory_Intf::GetW() {
  auto p = dynamic_cast<DataFactory_Impl<T, F>* >(this);
  if (!p) {
    Errors::Message msg;
    msg << "factory requested via incorrect type: \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->GetW();
}

} // namespace Amanzi

#endif
