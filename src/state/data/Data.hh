/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_DATA_HH_
#define AMANZI_DATA_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "Data_Impl.hh"

namespace Amanzi {

class Visualization;
class Checkpoint;

// Interface of a Data
//
// This interface allows Data objects to be kept in containers by
// providing all the necessary constructors, operator=, etc.  This is really
// a wrapper for std::unique_ptr -- equivalent code could directly store
// pointers to Data_Intf objects.
//
class Data {
 public:
  //
  // Default constructor
  //
  // This should never be used, only exists to make containers happy.  This is
  // NOT a valid object as it has no type information.
  Data() : p_(std::unique_ptr<Data_Intf>()) {}

  //
  // Constructor with type information
  //
  // This should not be used, instead use the non-member function:
  // CreateNullData<T>().
  Data(std::unique_ptr<Data_Intf> t) : p_(std::move(t)) {}

  // Copy constructor deleted, as we don't necessarily know how to copy
  // construct
  Data(const Data& other) = delete;

  // move constructor
  Data(Data&& other) noexcept : p_(std::move(other.p_)) {}

  // steal an r-value
  void swap(Data&& other) noexcept { p_.swap(other.p_); }

  // operator= with lvalue reference sets the values equal
  Data& operator=(const Data& other)
  {
    if (&other != this) {
      if (!p_) {
        Errors::Message msg;
        msg << " data not created through RecordSet::SetType() or "
               "State::CreatData()";
        throw(msg);
      }
      *p_ = *other.p_;
    }
    return *this;
  }

  // operator= with rvalue steals the data via swap
  Data& operator=(Data&& other) = default;

  // accessor -- const ref
  template <typename T>
  const T& Get() const
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    return p_->Get<T>();
  }

  // accessor -- non-const ref
  template <typename T>
  T& GetW()
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    return p_->GetW<T>();
  }

  // accessor -- const pointer
  template <typename T>
  Teuchos::RCP<const T> GetPtr() const
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    return p_->GetPtr<T>();
  }

  // accessor -- non-const shared pointer
  template <typename T>
  Teuchos::RCP<T> GetPtrW()
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    return p_->GetPtrW<T>();
  }

  // mutator -- set data by pointer
  template <typename T>
  void SetPtr(Teuchos::RCP<T> t)
  {
    if (!p_) { p_ = std::make_unique<Data_Impl<T>>(t); }
    p_->SetPtr(t);
  }

  // mutator -- set value
  template <typename T>
  void Set(const T& t)
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    p_->Set(t);
  }

  // virtual interface for ad-hoc polymorphism
  void
  WriteVis(const Visualization& vis, const Teuchos::ParameterList& attrs) const
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    p_->WriteVis(vis, attrs);
  }
  void WriteCheckpoint(const Checkpoint& chkp,
                       const Teuchos::ParameterList& attrs) const
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    p_->WriteCheckpoint(chkp, attrs);
  }
  void
  ReadCheckpoint(const Checkpoint& chkp, const Teuchos::ParameterList& attrs)
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    p_->ReadCheckpoint(chkp, attrs);
  }
  bool
  Initialize(Teuchos::ParameterList& plist, const Teuchos::ParameterList& attrs)
  {
    if (!p_) {
      Errors::Message msg;
      msg << " data not created through RecordSet::SetType() or "
             "State::CreatData()";
      throw(msg);
    }
    return p_->Initialize(plist, attrs);
  }

 private:
  std::unique_ptr<Data_Intf> p_;
};

//
// Non-member constructor of default (empty) Data
// --------------------------------------------------------
template <typename T>
Data
data()
{
  return Data(std::unique_ptr<Data_Intf>(new Data_Impl<T>()));
}

//
// Non-member constructor of Data with RCP
// --------------------------------------------------------
template <typename T>
Data
data(const Teuchos::RCP<T>& p)
{
  return Data(std::unique_ptr<Data_Intf>(new Data_Impl<T>(p)));
}

// //
// // Non-member constructor of Data with value type.
// // --------------------------------------------------------
// template<typename T, typename...Arg>
// Data data(Arg&&...arg) {
//   return Data(std::unique_ptr<Data_Intf>(new
//   Data_Impl<T>(std::forward<Arg...>(arg...))));
// }

} // namespace Amanzi

#endif
