/*
  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
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
namespace Impl {

template <typename T>
class HasCreate {
 private:
  HasCreate() = delete;

  struct one {
    char x[1];
  };
  struct two {
    char x[2];
  };

  template <typename C>
  static one test(decltype(void(std::declval<C&>().Create()))*);
  template <typename C>
  static two test(...);

 public:
  static constexpr bool value = sizeof(test<T>(0)) == sizeof(one);
};


// thing factory base class
class DataFactory_Intf {
 public:
  virtual ~DataFactory_Intf(){};

  virtual std::unique_ptr<DataFactory_Intf> Clone() const = 0;

  template <typename T, typename F>
  bool ValidType() const;

  virtual void Create(Data& t) const = 0;
  virtual Data Create() const = 0;

  template <typename T, typename F>
  const F& Get() const;

  template <typename T, typename F>
  F& GetW();
};


// thing factory implementation
template <typename T, typename F>
class DataFactory_Impl : public DataFactory_Intf {
 public:
  explicit DataFactory_Impl(const F& f) : f_(std::move(std::unique_ptr<F>(new F(f)))) {}

  DataFactory_Impl(const DataFactory_Impl& other)
    : f_(std::move(std::unique_ptr<F>(new F(*other.f_))))
  {}

  std::unique_ptr<DataFactory_Intf> Clone() const override
  {
    return std::unique_ptr<DataFactory_Impl<T, F>>(new DataFactory_Impl<T, F>(*this));
  }

  void Create(Data& t) const override
  {
    t = Create();
    return;
  }

  Data Create() const override { return Create_(); }

  const F& Get() const { return *f_; }

  F& GetW() { return *f_; }

 private:
  std::unique_ptr<F> f_;


  // two ways to Create a T from an F.  Either:
  // Constructor:  T(F);
  // Create method: Teuchos::RCP<T> = F.Create();
  //
  template <class Q = F>
  typename std::enable_if<HasCreate<Q>::value, Data>::type Create_() const
  {
    return data<T>(f_->Create());
  }

  template <class Q = F>
  typename std::enable_if<!HasCreate<Q>::value, Data>::type Create_() const
  {
    return data<T>(Teuchos::rcp(new T(*f_)));
  }
};


// partial specialization for NullFactory
template <typename T>
class DataFactory_Impl<T, NullFactory> : public DataFactory_Intf {
 public:
  explicit DataFactory_Impl(const NullFactory& f)
    : f_(std::unique_ptr<NullFactory>(new NullFactory(f))){};

  DataFactory_Impl(const DataFactory_Impl& other)
    : f_(std::unique_ptr<NullFactory>(new NullFactory(*other.f_))){};

  std::unique_ptr<DataFactory_Intf> Clone() const override
  {
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
template <typename T, typename F>
bool
DataFactory_Intf::ValidType() const
{
  auto p = dynamic_cast<const DataFactory_Impl<T, F>*>(this);
  if (!p) {
    return false;
  } else {
    return true;
  }
}

template <typename T, typename F>
const F&
DataFactory_Intf::Get() const
{
  auto p = dynamic_cast<DataFactory_Impl<T, F>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << "factory requested via incorrect type: \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->Get();
}

template <typename T, typename F>
F&
DataFactory_Intf::GetW()
{
  auto p = dynamic_cast<DataFactory_Impl<T, F>*>(this);
  if (!p) {
    Errors::Message msg;
    msg << "factory requested via incorrect type: \"" << typeid(T).name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return p->GetW();
}

} // namespace Impl
} // namespace Amanzi

#endif
