/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*

  CompVector

  Interface for CompVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.

*/

#ifndef AMANZI_COMPVECTOR_IMPL_HH_
#define AMANZI_COMPVECTOR_IMPL_HH_

#include <numeric>
#include "dbc.hh"
#include "errors.hh"

#include "AmanziMap.hh"
#include "AmanziVector.hh"

#include "VectorHarness.hh"

namespace Amanzi {

// Constructor
template<typename Scalar>
CompVector<Scalar>::CompVector(const Teuchos::RCP<const CompMap>& map)
    : map_(map) {};


// copy constructor
template<typename Scalar>
CompVector<Scalar>::CompVector(const CompVector<Scalar>& other)
    : CompVector(other.Map())
{
  CreateData();
};

// Check consistency of meta-data and allocate data.
template<typename Scalar>
void
CompVector<Scalar>::CreateData() {
#ifdef ENABLE_DBC
  if (data_[0] != Teuchos::null) {
    std::cout << "WARNING: CompositeVector::CreateData() called on already created vector!" << std::endl;
  }
#endif

  // create the data
  for (const auto& name : *this) {
    auto map = ComponentMap(name);
    if (map != Teuchos::null) {
      data_[name] = Teuchos::rcp(new MultiVector_type_<Scalar>(map, NumVectors(name), true));
    }
  }
};


// Assigment.
template<typename Scalar>
CompVector<Scalar>&
CompVector<Scalar>::operator=(const CompVector<Scalar>& other) {
  if (this != &other) {
    for (const auto& name : *this) {
      const auto other_v = other.GetComponent(name);
      const auto this_v = GetComponent(name);
      if (other_v == Teuchos::null) {
        SetComponent(Teuchos::null);
      } else if (this_v == Teuchos::null) {
        SetComponent(name, Tpetra::createCopy(*other_v));
      } else {
        Tpetra::deep_copy(*this_v, *other_v);
      }
    }
  }
  return *this;
};

template<typename Scalar>
LO
CompVector<Scalar>::MyLength(const std::string& name) const {
    return ComponentMap(name)->getNodeNumElements();
}


// View data, const version.
// I would prefer these be private, but for now...
template<typename Scalar>
cMultiVector_ptr_type_<Scalar>
CompVector<Scalar>::GetComponent(const std::string& name) const {
  if (!Map()->HasComponent(name)) {
    Errors::Message message("CompVector: Requested component ("+name+") does not exist.");
    throw(message);
  }
  return GetComponent_(name);
}


// View data, non-const version.
template<typename Scalar>
MultiVector_ptr_type_<Scalar>
CompVector<Scalar>::GetComponent(const std::string& name) {
  if (!Map()->HasComponent(name)) {
    Errors::Message message("CompVector: Requested component ("+name+") does not exist.");
    throw(message);
  }
  return GetComponent_(name);
};


// Set data
template<typename Scalar>
void
CompVector<Scalar>::SetComponent(const std::string& name,
        const MultiVector_ptr_type_<Scalar>& data)
{
  if (!Map()->HasComponent(name)) {
    Errors::Message message("CompVector: Requested SetComponent("+name+") does not exist.");
    throw(message);
  }
  if (data != Teuchos::null) {
    if (!ComponentMap(name)->isSameAs(*data->getMap())) {
      Errors::Message message("CompVector: Requested SetComponent("+name+") has incompatible map.");
      throw(message);
    }
    if (!(NumVectors(name) == data->getNumVectors())) {
      Errors::Message message("CompVector: Requested SetComponent("+name+") has incompatible number of DoFs.");
      throw(message);
    }
  }
  data_[name] = data;
}

// -- View a component vector
template<typename Scalar>
template<class DeviceType>
cMultiVectorView_type_<DeviceType,Scalar>
CompVector<Scalar>::ViewComponent(const std::string& name) const
{
  using memory_space = typename DeviceType::memory_space;
  return VectorHarness::getMultiVector(VectorHarness::readOnly(*GetComponent_(name)).on(memory_space()));
}

template<typename Scalar>
template<class DeviceType>
MultiVectorView_type_<DeviceType,Scalar>
CompVector<Scalar>::ViewComponent(const std::string& name)
{
  using memory_space = typename DeviceType::memory_space;
  return VectorHarness::getMultiVector(VectorHarness::readWrite(*GetComponent_(name)).on(memory_space()));
}

// -- SubView of a component vector
template<typename Scalar>
template<class DeviceType>
cVectorView_type_<DeviceType,Scalar>
CompVector<Scalar>::ViewComponent(const std::string& name, std::size_t dof) const
{
  using memory_space = typename DeviceType::memory_space;
  return Kokkos::subview(ViewComponent<DeviceType>(name), Kokkos::ALL(), dof);
}

template<typename Scalar>
template<class DeviceType>
VectorView_type_<DeviceType,Scalar>
CompVector<Scalar>::ViewComponent(const std::string& name, std::size_t dof)
{
  using memory_space = typename DeviceType::memory_space;
  return Kokkos::subview(ViewComponent<DeviceType>(name), Kokkos::ALL(), dof);
}



// Vector operations.
// -- Insert value into data.
template<typename Scalar>
int
CompVector<Scalar>::PutScalar(Scalar scalar) {
  int ierr = 0;
  for (const auto& name : *this) {
    GetComponent_(name)->putScalar(scalar);
  }
  return ierr;
};


// // -- Insert values into data, by DOF, not by component!
// template<typename Scalar>
// int
// CompVector<Scalar>::PutScalar(std::vector<Scalar> scalar) {
//   for (const auto& name : *this) {
//     AMANZI_ASSERT(scalar.size() == num_dofs_[i]);
//     for (int lcv_vector = 0; lcv_vector != GetComponent_(name)->getNumVectors(); ++lcv_vector) {
//       (*GetComponent_(name))(lcv_vector)->putScalar(scalar[lcv_vector]);
//     }
//   }
//   return 0;
// };


// -- Insert value into component [name].
template<typename Scalar>
int
CompVector<Scalar>::PutScalar(const std::string& name, Scalar scalar) {
  GetComponent(name)->putScalar(scalar);
  return 0;
};


// -- Insert values into data of component [name].
template<typename Scalar>
int
CompVector<Scalar>::PutScalar(const std::string& name, const std::vector<Scalar>& scalar) {
  if (scalar.size() != NumComponents(name)) {
    Errors::Message message("CompVector: requested PutScalar("+name+") with incorrectly sized vector of scalars.");
    throw(message);
  }

  for (int lcv_vector = 0; lcv_vector != NumVectors(name); ++lcv_vector) {
    (*GetComponent_(name))(lcv_vector)->putScalar(scalar[lcv_vector]);
  }
  return 0;
};


// this <- abs(this)
template<typename Scalar>
int
CompVector<Scalar>::Abs(const CompVector<Scalar>& other) {
  for (const auto& name : *this) {
    GetComponent_(name)->abs(*other.GetComponent_(name));
  }
  return 0;
}


// -- this <- value*this
template<typename Scalar>
int
CompVector<Scalar>::Scale(Scalar value) {
  for (const auto& name : *this) {
    GetComponent_(name)->scale(value);
  }
  return 0;
};


// Scale() applied to component name.
template<typename Scalar>
int
CompVector<Scalar>::Scale(const std::string& name, Scalar value) {
  GetComponent_(name)->scale(value);
  return 0;
};


// this <- abs(this)
template<typename Scalar>
int
CompVector<Scalar>::Reciprocal(const CompVector<Scalar>& other) {
  for (const auto& name : *this) {
    GetComponent_(name)->reciprocal(*other.GetComponent_(name));
  }
  return 0;
}


// -- result <- other \dot this
template<typename Scalar>
int
CompVector<Scalar>::Dot(const CompVector<Scalar>& other, Scalar* result) const {
  *result = 0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> intermediate_result(NumVectors(name));
    GetComponent_(name)->dot(*(other.GetComponent_(name)), intermediate_result);
    *result += std::accumulate(intermediate_result.begin(), intermediate_result.end(), 0);
  }
  return 0;
};


// -- this <- scalarA*A + scalarThis*this
template<typename Scalar>
CompVector<Scalar>&
CompVector<Scalar>::Update(Scalar scalarA, const CompVector<Scalar>& A, Scalar scalarThis) {
  for (const auto& name : *this) {
    GetComponent_(name)->update(scalarA, *A.GetComponent_(name), scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
template<typename Scalar>
CompVector<Scalar>&
CompVector<Scalar>::Update(Scalar scalarA, const CompVector<Scalar>& A,
                           Scalar scalarB, const CompVector<Scalar>& B, Scalar scalarThis) {
  for (const auto& name : *this) {
    GetComponent_(name)->update(scalarA, *A.GetComponent_(name), scalarB, *B.GetComponent_(name), scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
template<typename Scalar>
int
CompVector<Scalar>::Multiply(Scalar scalarAB, const CompVector<Scalar>& A, const CompVector<Scalar>& B,
                             Scalar scalarThis) {
  for (const auto& name : *this) {
    if (A.GetComponent_(name)->getNumVectors() > 1) {
      Errors::Message message("Not implemented multiply: Tpetra does not provide elementwise multiply.");
      Exceptions::amanzi_throw(message);
    }
    GetComponent_(name)->elementWiseMultiply(scalarAB, *A.GetComponent_(name)->getVector(0), *B.GetComponent_(name), scalarThis);
  }
  return 0;
};

// // -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise division
// template<typename Scalar>
// int
// CompVector<Scalar>::ReciprocalMultiply(Scalar scalarAB, const CompVector<Scalar>& A, const CompVector<Scalar>& B,
//                   Scalar scalarThis) {
//   for (const auto& name : *this) {
//     GetComponent_(name)->reciprocalMultiply(scalarAB, *A.GetComponent_(name), *B.GetComponent_(name), scalarThis);
//   }
//   return 0;
// };


// -- norms
template<typename Scalar>
int
CompVector<Scalar>::NormInf(Scalar* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  *norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->normInf(norm_locs());
    Scalar my_norm_loc = *std::max_element(norm_locs.begin(), norm_locs.end());
    *norm = std::max(my_norm_loc, *norm);
  }
  return 0;
};


template<typename Scalar>
int
CompVector<Scalar>::Norm1(Scalar* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  *norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->norm1(norm_locs());
    *norm += std::accumulate(norm_locs.begin(), norm_locs.end(), 0);
  }
  return 0;
};


template<typename Scalar>
int
CompVector<Scalar>::Norm2(Scalar* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  *norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->norm2(norm_locs());
    *norm += std::accumulate(norm_locs.begin(), norm_locs.end(), 0, [](Scalar x, Scalar y) { return x + y*y; } );
  }
  *norm = sqrt(*norm);
  return 0;
};



// template<typename Scalar>
// int
// CompVector<Scalar>::MinValue(Scalar* value) const {
//   if (value == NULL) return 1;
//   if (data_.size() == 0) return 1;

//   int ierr = 0;
//   Scalar value_loc[1];

//   *value = 1e+50;
//   for (const auto& name : *this) {
//     for (int lcv_vector = 0; lcv_vector != GetComponent_(name)->NumVectors(); ++lcv_vector) {
//       ierr = (*GetComponent_(name))(lcv_vector)->MinValue(value_loc);
//       if (ierr) return ierr;
//       *value = std::min(*value, value_loc[0]);
//     }
//   }
//   return ierr;
// };


// template<typename Scalar>
// int
// CompVector<Scalar>::MaxValue(Scalar* value) const {
//   if (value == NULL) return 1;
//   if (data_.size() == 0) return 1;

//   int ierr = 0;
//   Scalar value_loc[1];

//   *value = -1e+50;
//   for (const auto& name : *this) {
//     for (int lcv_vector = 0; lcv_vector != GetComponent_(name)->NumVectors(); ++lcv_vector) {
//       ierr = (*GetComponent_(name))(lcv_vector)->MaxValue(value_loc);
//       if (ierr) return ierr;
//       *value = std::max(*value, value_loc[0]);
//     }
//   }
//   return ierr;
// };


// template<typename Scalar>
// int
// CompVector<Scalar>::MeanValue(Scalar* value) const {
//   if (value == NULL) return 1;
//   if (data_.size() == 0) return 1;

//   int ierr(0), n(0), n_loc;
//   Scalar value_loc[1];

//   *value = 0.0;
//   for (const auto& name : *this) {
//     n_loc = GetComponent_(name)->GlobalLength(); 
//     for (int lcv_vector = 0; lcv_vector != GetComponent_(name)->NumVectors(); ++lcv_vector) {
//       ierr = (*GetComponent_(name))(lcv_vector)->MeanValue(value_loc);
//       if (ierr) return ierr;
//       *value += value_loc[0] * n_loc;
//       n += n_loc;
//     }
//   }
//   *value /= n;
//   return ierr;
// };


// Debugging?
template<typename Scalar>
void
CompVector<Scalar>::Print(std::ostream& os, bool data_io) const {
  os << "Comp Vector" << std::endl;
  os << "  components: ";
  for (const auto& name : *this) {
    os << name << "(" << GetComponent_(name)->getNumVectors() << ") ";
  }
  os << std::endl;
  if (data_io) {
    for (const auto& name : *this) {
      GetComponent_(name)->print(os);
    }
  }
};


// Populate by random numbers between -1 and 1. 
template<typename Scalar>
int
CompVector<Scalar>::Random() {
  for (const auto& name : *this) {
    GetComponent_(name)->randomize();
  }
  return 0;
};


} // namespace

#endif

