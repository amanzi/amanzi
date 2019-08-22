/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*

  BlockVector

  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.

*/

#ifndef AMANZI_BLOCK_VECTOR_IMPL_HH_
#define AMANZI_BLOCK_VECTOR_IMPL_HH_

#include <numeric>
#include "dbc.hh"
#include "errors.hh"

#include "AmanziMap.hh"
#include "AmanziVector.hh"

#include "VectorHarness.hh"

namespace Amanzi {

// Constructor
template<typename Scalar>
BlockVector<Scalar>::BlockVector(const Teuchos::RCP<const BlockSpace>& map, InitMode mode)
    : map_(map) {
  for (const auto& name : *this) {
    SetComponent_(name, true, Teuchos::null);
    SetComponent_(name, false, Teuchos::null);
  }
  CreateData_(mode);
};


// copy constructor
template<typename Scalar>
BlockVector<Scalar>::BlockVector(const BlockVector<Scalar>& other, InitMode mode)
    : BlockVector(other.Map())
{
  for (const auto& name : *this) {
    SetComponent_(name, true, Teuchos::null);
    SetComponent_(name, false, Teuchos::null);
  }
  CreateData_(mode);
  if (mode == InitMode::COPY) {
    *this = other;
  }
};

// Check consistency of meta-data and allocate data.
template<typename Scalar>
void
BlockVector<Scalar>::CreateData_(InitMode mode) {
  // NOTE: Tpetra zero-initializes by default, so there is no difference between NONE and ZERO
  if (mode != InitMode::NOALLOC) {
    // create the data
    for (const auto& name : *this) {
      auto ghost_map = Map()->ComponentMap(name, true);
      if (ghost_map != Teuchos::null) {
        SetComponent_(name, true, Teuchos::rcp(new MultiVector_type_<Scalar>(ghost_map, NumVectors(name), true)));

        auto g_comp = GetComponent_(name, true);
        auto master_map = Map()->ComponentMap(name, false);

        // get the ghost component's data
        auto m_comp = g_comp->offsetViewNonConst(master_map, 0);
        SetComponent_(name, false, m_comp);

      } else {
        auto master_map = Map()->ComponentMap(name, false);
        SetComponent_(name, false, Teuchos::rcp(new MultiVector_type_<Scalar>(master_map, NumVectors(name), true)));
      }    
    }
  }
};


// Assigment.
template<typename Scalar>
BlockVector<Scalar>&
BlockVector<Scalar>::operator=(const BlockVector<Scalar>& other) {
  if (this != &other) {
    for (const auto& name : *this) {
      const auto other_v = other.GetComponent(name, false);
      const auto this_v = GetComponent(name, false);
      if (other_v == Teuchos::null && this_v == Teuchos::null) {
        // pass
      } else if (this_v == Teuchos::null) {
        auto new_v = Teuchos::rcp(new MultiVector_type_<Scalar>(*other_v));
        SetComponent_(name, false, new_v);
      } else {
        Tpetra::deep_copy(*this_v, *other_v);
      }
    }
  }
  return *this;
};


// View data, const version.
// I would prefer these be private, but for now...
template<typename Scalar>
cMultiVector_ptr_type_<Scalar>
BlockVector<Scalar>::GetComponent(const std::string& name, bool ghosted) const {
  if (!Map()->HasComponent(name)) {
    Errors::Message message("BlockVector: Requested component ("+name+") does not exist.");
    throw(message);
  }
  return GetComponent_(name, ghosted);
}


// View data, non-const version.
template<typename Scalar>
MultiVector_ptr_type_<Scalar>
BlockVector<Scalar>::GetComponent(const std::string& name, bool ghosted) {
  if (!Map()->HasComponent(name)) {
    Errors::Message message("BlockVector: Requested component ("+name+") does not exist.");
    throw(message);
  }
  return GetComponent_(name, ghosted);
};


// // Set data
// template<typename Scalar>
// void
// BlockVector<Scalar>::SetComponent(const std::string& name,
//         const MultiVector_ptr_type_<Scalar>& data)
// {
//   if (!Map()->HasComponent(name)) {
//     Errors::Message message("BlockVector: Requested SetComponent("+name+") does not exist.");
//     throw(message);
//   }
//   if (data != Teuchos::null) {
//     if (!ComponentMap(name)->isSameAs(*data->getMap())) {
//       Errors::Message message("BlockVector: Requested SetComponent("+name+") has incompatible map.");
//       throw(message);
//     }
//     if (!(NumVectors(name) == data->getNumVectors())) {
//       Errors::Message message("BlockVector: Requested SetComponent("+name+") has incompatible number of DoFs.");
//       throw(message);
//     }
//   }
//   data_[name] = data;
// }

// -- View a component vector
template<typename Scalar>
template<class DeviceType>
cMultiVectorView_type_<DeviceType,Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, bool ghosted) const
{
  using memory_space = typename DeviceType::memory_space;
  return VectorHarness::getMultiVector(VectorHarness::readOnly(*GetComponent_(name, ghosted)).on(memory_space()));
}

template<typename Scalar>
template<class DeviceType>
MultiVectorView_type_<DeviceType,Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, bool ghosted)
{
  using memory_space = typename DeviceType::memory_space;
  return VectorHarness::getMultiVector(VectorHarness::readWrite(*GetComponent_(name, ghosted)).on(memory_space()));
}

// -- SubView of a component vector
template<typename Scalar>
template<class DeviceType>
cVectorView_type_<DeviceType,Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, std::size_t dof, bool ghosted) const
{
  using memory_space = typename DeviceType::memory_space;
  return Kokkos::subview(ViewComponent<DeviceType>(name, ghosted), Kokkos::ALL(), dof);
}

template<typename Scalar>
template<class DeviceType>
VectorView_type_<DeviceType,Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, std::size_t dof, bool ghosted)
{
  using memory_space = typename DeviceType::memory_space;
  return Kokkos::subview(ViewComponent<DeviceType>(name, ghosted), Kokkos::ALL(), dof);
}


// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
template<typename Scalar>
void
BlockVector<Scalar>::ScatterMasterToGhosted(Tpetra::CombineMode mode) const
{
  for (const auto& name : *this) {
    ScatterMasterToGhosted(name, mode);
  }
};


// Scatter master values to ghosted values, on all components, in a mode.
//
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
//
// Note that although scatter changes things, it doesn't change master
// data, so we allow it to work on const.  This is necessary for a
// non-owning PK to communicate a non-owned vector.
//
// This Scatter() is not managed, and is always done.  Tags changed.
template<typename Scalar>
void
BlockVector<Scalar>::ScatterMasterToGhosted(const std::string& name, Tpetra::CombineMode mode) const
{
  auto import = Map()->Importer(name);
  if (import == Teuchos::null) return;
  auto g_comp = ghost_data_.at(name); // note cannot use GetComponent_() here, need non-const
  if (g_comp == Teuchos::null) return;
  auto m_comp = GetComponent_(name, false);
  g_comp->doImport(*m_comp, *import, mode);
};




// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
template<typename Scalar>
void
BlockVector<Scalar>::GatherGhostedToMaster(Tpetra::CombineMode mode)
{
  for (const auto& name : *this) {
    GatherGhostedToMaster(name, mode);
  }
};


template<typename Scalar>
void
BlockVector<Scalar>::GatherGhostedToMaster(const std::string& name, Tpetra::CombineMode mode)
{
  auto import = Map()->Importer(name);
  if (import == Teuchos::null) return;
  // communicate
  auto g_comp = GetComponent(name, true);
  auto m_comp = GetComponent(name, false);
  m_comp->doExport(*g_comp, *import, mode);
};


// Vector operations.
// -- Insert value into data.
template<typename Scalar>
int
BlockVector<Scalar>::PutScalar(Scalar scalar) {
  int ierr = 0;
  for (const auto& name : *this) {
    GetComponent_(name)->putScalar(scalar);
  }
  return ierr;
};


// -- Insert value into component [name].
template<typename Scalar>
int
BlockVector<Scalar>::PutScalar(const std::string& name, Scalar scalar) {
  GetComponent(name)->putScalar(scalar);
  return 0;
};


// -- Insert values into data of component [name].
template<typename Scalar>
int
BlockVector<Scalar>::PutScalar(const std::string& name, const std::vector<Scalar>& scalar) {
  if (scalar.size() != NumVectors(name)) {
    Errors::Message message("BlockVector: requested PutScalar("+name+") with incorrectly sized vector of scalars.");
    throw(message);
  }

  for (int lcv_vector = 0; lcv_vector != NumVectors(name); ++lcv_vector) {
    (*GetComponent_(name))(lcv_vector)->putScalar(scalar[lcv_vector]);
  }
  return 0;
};


// Vector operations.
// -- Insert value into data.
template<typename Scalar>
int
BlockVector<Scalar>::PutScalarMasterAndGhosted(Scalar scalar) {
  int ierr = 0;
  for (const auto& name : *this) {
    GetComponent_(name,true)->putScalar(scalar);
  }
  return ierr;
};


// template<typename Scalar>
// int
// BlockVector<Scalar>::PutScalarGhosted(Scalar scalar) {
//   for (int lcv_comp = 0; lcv_comp != NumComponents(); ++lcv_comp) {
//     int size_owned = mastervec_->size(names_[lcv_comp]);
//     int size_ghosted = ghostvec_->size(names_[lcv_comp]);

//     using Range_type = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, int>;
//     /*
//     auto vec = ViewComponent(names_[lcv_comp], true);
//     Kokkos::parallel_for("PutScalarGhosted", Range_type(size_owned, size_ghosted), 
// 			 [=] (const int& i) { vec(i) = scalar; })
//     */
//     auto vec = ViewComponent(names_[lcv_comp], true);
//     auto vec_ghost_view = Kokkos::subview(vec, std::make_pair(size_owned, size_ghosted), Kokkos::ALL());
//     Kokkos::parallel_for("PutScalarGhosted", Range_type(0, size_ghosted - size_owned), 
    

//   }
//   return 0;
// }


// this <- abs(this)
template<typename Scalar>
int
BlockVector<Scalar>::Abs(const BlockVector<Scalar>& other) {
  for (const auto& name : *this) {
    GetComponent_(name)->abs(*other.GetComponent_(name));
  }
  return 0;
}


// -- this <- value*this
template<typename Scalar>
int
BlockVector<Scalar>::Scale(Scalar value) {
  for (const auto& name : *this) {
    GetComponent_(name)->scale(value);
  }
  return 0;
};


// Scale() applied to component name.
template<typename Scalar>
int
BlockVector<Scalar>::Scale(const std::string& name, Scalar value) {
  GetComponent_(name)->scale(value);
  return 0;
};


// this <- abs(this)
template<typename Scalar>
int
BlockVector<Scalar>::Reciprocal(const BlockVector<Scalar>& other) {
  for (const auto& name : *this) {
    GetComponent_(name)->reciprocal(*other.GetComponent_(name));
  }
  return 0;
}


// -- result <- other \dot this
template<typename Scalar>
int
BlockVector<Scalar>::Dot(const BlockVector<Scalar>& other, Scalar* result) const {
  *result = 0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> intermediate_result(NumVectors(name));
    GetComponent_(name)->dot(*(other.GetComponent_(name)), intermediate_result);
    *result += std::accumulate(intermediate_result.begin(), intermediate_result.end(), (Scalar) 0);
  }
  return 0;
};


// -- this <- scalarA*A + scalarThis*this
template<typename Scalar>
BlockVector<Scalar>&
BlockVector<Scalar>::Update(Scalar scalarA, const BlockVector<Scalar>& A, Scalar scalarThis) {
  for (const auto& name : *this) {
    GetComponent_(name)->update(scalarA, *A.GetComponent_(name), scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
template<typename Scalar>
BlockVector<Scalar>&
BlockVector<Scalar>::Update(Scalar scalarA, const BlockVector<Scalar>& A,
                           Scalar scalarB, const BlockVector<Scalar>& B, Scalar scalarThis) {
  for (const auto& name : *this) {
    GetComponent_(name)->update(scalarA, *A.GetComponent_(name), scalarB, *B.GetComponent_(name), scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
template<typename Scalar>
int
BlockVector<Scalar>::Multiply(Scalar scalarAB, const BlockVector<Scalar>& A, const BlockVector<Scalar>& B,
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
// BlockVector<Scalar>::ReciprocalMultiply(Scalar scalarAB, const BlockVector<Scalar>& A, const BlockVector<Scalar>& B,
//                   Scalar scalarThis) {
//   for (const auto& name : *this) {
//     GetComponent_(name)->reciprocalMultiply(scalarAB, *A.GetComponent_(name), *B.GetComponent_(name), scalarThis);
//   }
//   return 0;
// };


// -- norms
template<typename Scalar>
int
BlockVector<Scalar>::NormInf(Scalar* norm) const {
  if (norm == NULL) return 1;
  if (master_data_.size() == 0) return 1;

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
BlockVector<Scalar>::Norm1(Scalar* norm) const {
  if (norm == NULL) return 1;
  if (master_data_.size() == 0) return 1;

  *norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->norm1(norm_locs());
    *norm += std::accumulate(norm_locs.begin(), norm_locs.end(), (Scalar) 0);
  }
  return 0;
};


template<typename Scalar>
int
BlockVector<Scalar>::Norm2(Scalar* norm) const {
  if (norm == NULL) return 1;
  if (master_data_.size() == 0) return 1;

  *norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->norm2(norm_locs());
    *norm += std::accumulate(norm_locs.begin(), norm_locs.end(), (Scalar) 0, [](Scalar x, Scalar y) { return x + y*y; } );
  }
  *norm = sqrt(*norm);
  return 0;
};



// template<typename Scalar>
// int
// BlockVector<Scalar>::MinValue(Scalar* value) const {
//   if (value == NULL) return 1;
//   if (master_data_.size() == 0) return 1;

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
// BlockVector<Scalar>::MaxValue(Scalar* value) const {
//   if (value == NULL) return 1;
//   if (master_data_.size() == 0) return 1;

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
// BlockVector<Scalar>::MeanValue(Scalar* value) const {
//   if (value == NULL) return 1;
//   if (master_data_.size() == 0) return 1;

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
BlockVector<Scalar>::Print(std::ostream& os, bool ghosted, bool data_io) const {
  os << "Comp Vector" << std::endl;
  os << "  components: ";
  for (const auto& name : *this) {
    os << name << "(" << NumVectors(name) << ") ";
  }
  os << std::endl;
  if (data_io) {
    for (const auto& name : *this) {
      GetComponent_(name, ghosted)->print(os);
    }
  }
};


// Populate by random numbers between -1 and 1. 
template<typename Scalar>
int
BlockVector<Scalar>::Random() {
  for (const auto& name : *this) {
    GetComponent_(name)->randomize();
  }
  return 0;
};


} // namespace

#endif

