/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

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

// thing factory wrapper
class DataFactory {
 public:
  DataFactory() : p_(std::unique_ptr<DataFactory_Intf>()) {};

  DataFactory(DataFactory_Intf *t) : p_(t) {}

  DataFactory(const DataFactory& other) : p_(other.p_->Clone()) {}

  DataFactory(DataFactory&& other) noexcept : p_(std::move(other.p_)) {}

  void swap(DataFactory& other) noexcept { p_.swap(other.p_); }

  DataFactory& operator=(DataFactory other) {
    if (&other != this)
      swap(other);
    return *this;
  }

  bool HasType() const { return p_.get(); }

  template <typename T, typename F> bool ValidType() const {
    return p_->ValidType<T, F>();
  }

  template <typename T, typename F> const F& Get() const {
    return p_->Get<T, F>();
  }

  template <typename T, typename F> F& GetW() { return p_->GetW<T, F>(); }

  void Create(Data& t) { return p_->Create(t); }

  Data Create() { return p_->Create(); }

 private:
  std::unique_ptr<DataFactory_Intf> p_;
};


template <typename T, typename F> DataFactory dataFactory() {
  return DataFactory(new DataFactory_Impl<T, F>());
}

template <typename T, typename F> DataFactory dataFactory(F f) {
  return DataFactory(new DataFactory_Impl<T, F>(std::move(f)));
}

} // namespace Amanzi

#endif
