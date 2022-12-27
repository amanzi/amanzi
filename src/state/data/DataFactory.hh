/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Basic factory that may store heterogeneous data and create
  heterogeneous Data objects to put into State.

  This shares the same design pattern as Data to do type erasure.
*/

#ifndef AMANZI_STATE_STOREABLE_FACTORY_HH_
#define AMANZI_STATE_STOREABLE_FACTORY_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "Data.hh"
#include "StateDefs.hh"

#include "DataFactory_Impl.hh"

namespace Amanzi {
namespace Impl {

// thing factory wrapper
class DataFactory {
 public:
  DataFactory() : p_(std::unique_ptr<DataFactory_Intf>()){};

  DataFactory(DataFactory_Intf* t, Data_Intf* d) : p_(t), data_p_(d) {}

  DataFactory(const DataFactory& other)
    : p_(other.p_->Clone()), data_p_(other.data_p_->CloneEmpty())
  {}

  DataFactory(DataFactory&& other) noexcept : p_(std::move(other.p_)) {}

  void swap(DataFactory& other) noexcept
  {
    p_.swap(other.p_);
    data_p_.swap(other.data_p_);
  }

  DataFactory& operator=(DataFactory other)
  {
    if (&other != this) swap(other);
    return *this;
  }

  bool HasType() const { return p_.get(); }

  template <typename T>
  bool ValidType() const
  {
    return data_p_->ValidType<T>();
  }

  template <typename T, typename F>
  bool ValidType() const
  {
    return p_->ValidType<T, F>();
  }

  template <typename T, typename F>
  const F& Get() const
  {
    return p_->Get<T, F>();
  }

  template <typename T, typename F>
  F& GetW()
  {
    return p_->GetW<T, F>();
  }

  void Create(Data& t) { return p_->Create(t); }

  Data Create() { return p_->Create(); }

 private:
  std::unique_ptr<DataFactory_Intf> p_;

  // note, this only exists to allow checking of the data type without knowing
  // the factory type.  It is not used outside of that, and data_p->p_ is
  // always null (it stores no data, just the type info).
  std::unique_ptr<Data_Intf> data_p_;
};


template <typename T, typename F>
DataFactory
dataFactory(const F& f)
{
  return DataFactory(new DataFactory_Impl<T, F>(f), new Data_Impl<T>());
}

template <typename T, typename F>
typename std::enable_if<std::is_default_constructible<F>::value, DataFactory>::type
dataFactory()
{
  F f;
  return DataFactory(new DataFactory_Impl<T, F>(f), new Data_Impl<T>());
}


// template <typename T, typename F> DataFactory dataFactory(F f) {
//   return DataFactory(new DataFactory_Impl<T, F>(std::move(f)));
// }

} // namespace Impl
} // namespace Amanzi

#endif
