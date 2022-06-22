/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

/*

  BlockVector

  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.

*/

#ifndef AMANZI_BLOCK_VECTOR_IMPL_HH_
#define AMANZI_BLOCK_VECTOR_IMPL_HH_

#include <numeric>
#include "Teuchos_FancyOStream.hpp"
#include "dbc.hh"
#include "errors.hh"

#include "AmanziMap.hh"
#include "AmanziVector.hh"


namespace Amanzi {

// Constructor
template <typename Scalar>
BlockVector<Scalar>::BlockVector(const Teuchos::RCP<const BlockSpace>& map,
                                 InitMode mode)
  : map_(map)
{
  for (const auto& name : *this) {
    SetComponent_(name, true, Teuchos::null);
    SetComponent_(name, false, Teuchos::null);
  }
  CreateData_(mode);
};


// copy constructor
template <typename Scalar>
BlockVector<Scalar>::BlockVector(const BlockVector<Scalar>& other,
                                 Teuchos::DataAccess access, InitMode mode)
  : map_(other.getMap())
{
  for (const auto& name : *this) {
    SetComponent_(name, true, Teuchos::null);
    SetComponent_(name, false, Teuchos::null);
  }

  if (access == Teuchos::DataAccess::View) {
    Errors::Message message("BlockVector: View semantic not supported.");
    throw(message);

    // for (const auto& name : *this) {
    //   SetComponent_(name, true, other.GetComponent_(name, true));
    //   SetComponent_(name, false, other.GetComponent_(name, false));
    // }
  } else {
    CreateData_(mode);
    if (mode == InitMode::COPY) { 
      *this = other; 
    }
  }
};

// Check consistency of meta-data and allocate data.
template <typename Scalar>
void
BlockVector<Scalar>::CreateData_(InitMode mode)
{
  // NOTE: Tpetra zero-initializes by default, so there is no difference between
  // NONE and ZERO
  if (mode != InitMode::NOALLOC) {
    // create the data
    for (const auto& name : *this) {
      auto ghost_map = getMap()->ComponentMap(name, true);
      if (ghost_map != Teuchos::null) {
        SetComponent_(name,
                      true,
                      Teuchos::rcp(new MultiVector_type_<Scalar>(
                        ghost_map, getNumVectors(name), true)));

        auto g_comp = GetComponent_(name, true);
        auto master_map = getMap()->ComponentMap(name, false);

        // get the ghost component's data
        auto m_comp = g_comp->offsetViewNonConst(master_map, 0);
        SetComponent_(name, false, m_comp);

      } else {
        auto master_map = getMap()->ComponentMap(name, false);
        SetComponent_(name,
                      false,
                      Teuchos::rcp(new MultiVector_type_<Scalar>(
                        master_map, getNumVectors(name), true)));
      }
    }
  }
};


// Assigment.
template <typename Scalar>
BlockVector<Scalar>&
BlockVector<Scalar>::operator=(const BlockVector<Scalar>& other)
{
  if (this != &other) {
    for (const auto& name : *this) {
      const auto other_v = other.GetComponent(name, false);
      const auto this_v = GetComponent(name, false);
      if (other_v == Teuchos::null && this_v == Teuchos::null) {
        // pass
      } else if (this_v == Teuchos::null) {
        auto new_v = Teuchos::rcp(new MultiVector_type_<Scalar>(*other_v,Teuchos::Copy));
        Tpetra::deep_copy(*new_v, *other_v);
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
template <typename Scalar>
cMultiVector_ptr_type_<Scalar>
BlockVector<Scalar>::GetComponent(const std::string& name, bool ghosted) const
{
  if (!getMap()->HasComponent(name)) {
    Errors::Message message("BlockVector: Requested component (" + name +
                            ") does not exist.");
    throw(message);
  }
  return GetComponent_(name, ghosted);
}


// View data, non-const version.
template <typename Scalar>
MultiVector_ptr_type_<Scalar>
BlockVector<Scalar>::GetComponent(const std::string& name, bool ghosted)
{
  if (!getMap()->HasComponent(name)) {
    Errors::Message message("BlockVector: Requested component (" + name +
                            ") does not exist.");
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
//   if (!getMap()->HasComponent(name)) {
//     Errors::Message message("BlockVector: Requested SetComponent("+name+")
//     does not exist."); throw(message);
//   }
//   if (data != Teuchos::null) {
//     if (!ComponentMap(name)->isSameAs(*data->getMap())) {
//       Errors::Message message("BlockVector: Requested SetComponent("+name+")
//       has incompatible map."); throw(message);
//     }
//     if (!(getNumVectors(name) == data->getNumVectors())) {
//       Errors::Message message("BlockVector: Requested SetComponent("+name+")
//       has incompatible number of DoFs."); throw(message);
//     }
//   }
//   data_[name] = data;
// }

// -- View a component vector
template <typename Scalar>
template <class DeviceType>
cMultiVectorView_type_<DeviceType, Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, bool ghosted) const
{
  using memory_space = typename DeviceType::memory_space;
  return GetComponent_(name,ghosted)->template getLocalView<memory_space>(Tpetra::Access::ReadOnly);
}

template <typename Scalar>
template <class DeviceType>
MultiVectorView_type_<DeviceType, Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, bool ghosted)
{
  using memory_space = typename DeviceType::memory_space;
  return  GetComponent_(name,ghosted)->template getLocalView<memory_space>(Tpetra::Access::ReadWrite);
}

// -- SubView of a component vector
template <typename Scalar>
template <class DeviceType>
cVectorView_type_<DeviceType, Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, std::size_t dof,
                                   bool ghosted) const
{
  return Kokkos::subview(
    ViewComponent<DeviceType>(name, ghosted), Kokkos::ALL(), dof);
}

template <typename Scalar>
template <class DeviceType>
VectorView_type_<DeviceType, Scalar>
BlockVector<Scalar>::ViewComponent(const std::string& name, std::size_t dof,
                                   bool ghosted)
{
  return Kokkos::subview(
    ViewComponent<DeviceType>(name, ghosted), Kokkos::ALL(), dof);
}


// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
template <typename Scalar>
void
BlockVector<Scalar>::ScatterMasterToGhosted(Tpetra::CombineMode mode) const
{
  for (const auto& name : *this) { ScatterMasterToGhosted(name, mode); }
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
template <typename Scalar>
void
BlockVector<Scalar>::ScatterMasterToGhosted(const std::string& name,
                                            Tpetra::CombineMode mode) const
{
  auto import = getMap()->Importer(name);
  if (import == Teuchos::null) return;
  auto g_comp = ghost_data_.at(
    name); // note cannot use GetComponent_() here, need non-const
  if (g_comp == Teuchos::null) return;
  auto m_comp = GetComponent_(name, false);
  g_comp->doImport(*m_comp, *import, mode);
};


// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
template <typename Scalar>
void
BlockVector<Scalar>::GatherGhostedToMaster(Tpetra::CombineMode mode)
{
  for (const auto& name : *this) { GatherGhostedToMaster(name, mode); }
};


template <typename Scalar>
void
BlockVector<Scalar>::GatherGhostedToMaster(const std::string& name,
                                           Tpetra::CombineMode mode)
{
  auto import = getMap()->Importer(name);
  if (import == Teuchos::null) return;
  // communicate
  auto g_comp = GetComponent(name, true);
  auto m_comp = GetComponent(name, false);
  m_comp->doExport(*g_comp, *import, mode);
};


// Vector operations.
// -- Insert value into data.
template <typename Scalar>
void
BlockVector<Scalar>::putScalar(Scalar scalar)
{
  for (const auto& name : *this) { GetComponent_(name)->putScalar(scalar); }
};


// -- Insert value into component [name].
template <typename Scalar>
void
BlockVector<Scalar>::putScalar(const std::string& name, Scalar scalar)
{
  GetComponent(name)->putScalar(scalar);
};


// -- Insert values into data of component [name].
template <typename Scalar>
void
BlockVector<Scalar>::putScalar(const std::string& name,
                               const std::vector<Scalar>& scalar)
{
  if (scalar.size() != getNumVectors(name)) {
    Errors::Message message("BlockVector: requested putScalar(" + name +
                            ") with incorrectly sized vector of scalars.");
    throw(message);
  }

  for (int lcv_vector = 0; lcv_vector != getNumVectors(name); ++lcv_vector) {
    (*GetComponent_(name))(lcv_vector)->putScalar(scalar[lcv_vector]);
  }
};


// Vector operations.
// -- Insert value into data.
template <typename Scalar>
void
BlockVector<Scalar>::putScalarMasterAndGhosted(Scalar scalar)
{
  for (const auto& name : *this) {
    GetComponent_(name, true)->putScalar(scalar);
  }
};


template<typename Scalar>
void
BlockVector<Scalar>::putScalarGhosted(Scalar scalar) {
  for (const auto& comp : *this) {
    auto vv = ViewComponent(comp, true);

    int size_owned = GetComponent(comp,false)->getLocalLength();
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> range_policy(
        {size_owned,0}, {vv.extent(0),vv.extent(1)});
    Kokkos::parallel_for(
        "BlockVector_impl::putScalarGhosted",
        range_policy,
        KOKKOS_LAMBDA (const int i, const int j) {
          vv(i,j) = scalar;
        });
  }
}


// this <- abs(this)
template <typename Scalar>
void
BlockVector<Scalar>::abs(const BlockVector<Scalar>& other)
{
  for (const auto& name : *this) {
    GetComponent_(name)->abs(*other.GetComponent_(name));
  }
}


// -- this <- value*this
template <typename Scalar>
void
BlockVector<Scalar>::scale(Scalar value)
{
  for (const auto& name : *this) { GetComponent_(name)->scale(value); }
};


// scale() applied to component name.
template <typename Scalar>
void
BlockVector<Scalar>::scale(const std::string& name, Scalar value)
{
  GetComponent_(name)->scale(value);
};


// this <- abs(this)
template <typename Scalar>
void
BlockVector<Scalar>::reciprocal(const BlockVector<Scalar>& other)
{
  for (const auto& name : *this) {
    GetComponent_(name)->reciprocal(*other.GetComponent_(name));
  }
}


// -- result <- other \dot this
template <typename Scalar>
Scalar
BlockVector<Scalar>::dot(const BlockVector<Scalar>& other) const
{
  Scalar result = 0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> intermediate_result(getNumVectors(name));
    GetComponent_(name)->dot(*(other.GetComponent_(name)), intermediate_result);
    result += std::accumulate(
      intermediate_result.begin(), intermediate_result.end(), (Scalar)0);
  }
  return result;
};


// -- this <- scalarA*A + scalarThis*this
template <typename Scalar>
void
BlockVector<Scalar>::update(Scalar scalarA, const BlockVector<Scalar>& A,
                            Scalar scalarThis)
{
  for (const auto& name : *this) {
    if (A.HasComponent(name)) {
      GetComponent_(name)->update(scalarA, *A.GetComponent_(name), scalarThis);
    }
  }
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
template <typename Scalar>
void
BlockVector<Scalar>::update(Scalar scalarA, const BlockVector<Scalar>& A,
                            Scalar scalarB, const BlockVector<Scalar>& B,
                            Scalar scalarThis)
{
  for (const auto& name : *this) {
    if (A.HasComponent(name) && B.HasComponent(name)) {
      GetComponent_(name)->update(scalarA,
              *A.GetComponent_(name),
              scalarB,
              *B.GetComponent_(name),
              scalarThis);
    }
  }
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
template <typename Scalar>
void
BlockVector<Scalar>::elementWiseMultiply(Scalar scalarAB,
                                         const BlockVector<Scalar>& A,
                                         const BlockVector<Scalar>& B,
                                         Scalar scalarThis)
{
  for (const auto& name : *this) {

    // There is a nasty gotcha here -- if this is A, when we call getVector we
    // clobber data and break things without erroring.  But this being B is ok.
    //
    if (B.GetComponent_(name) != this->GetComponent(name) &&
        A.GetComponent_(name) == this->GetComponent(name)) {
      if (B.GetComponent_(name)->getNumVectors() > 1) {
        Errors::Message message("Not implemented multiply: Tpetra does not "
                "provide multivector elementwise multiply.");
        Exceptions::amanzi_throw(message);
      }
    
      GetComponent_(name)->elementWiseMultiply(
          scalarAB,
          *B.GetComponent_(name)->getVector(0),
          *A.GetComponent_(name),
          scalarThis);
    } else if (B.GetComponent_(name) == this->GetComponent(name) &&
               A.GetComponent_(name) == this->GetComponent(name)) {
      Errors::Message message("Unclear whether, in elementWiseMultiply(), A,B,& this can all be the same vector.");
      Exceptions::amanzi_throw(message);
    } else {
      if (A.GetComponent_(name)->getNumVectors() > 1) {
        Errors::Message message("Not implemented multiply: Tpetra does not "
                "provide multivector elementwise multiply.");
        Exceptions::amanzi_throw(message);
      }
    
      GetComponent_(name)->elementWiseMultiply(
          scalarAB,
          *A.GetComponent_(name)->getVector(0),
          *B.GetComponent_(name),
          scalarThis);
    }
  }
};

// // -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise
// division template<typename Scalar> int
// BlockVector<Scalar>::ReciprocalelementWiseMultiply(Scalar scalarAB, const
// BlockVector<Scalar>& A, const BlockVector<Scalar>& B,
//                   Scalar scalarThis) {
//   for (const auto& name : *this) {
//     GetComponent_(name)->reciprocalelementWiseMultiply(scalarAB,
//     *A.GetComponent_(name), *B.GetComponent_(name), scalarThis);
//   }
//   return 0;
// };


// -- norms
template <typename Scalar>
Scalar
BlockVector<Scalar>::normInf() const
{
  Scalar norm = 0.;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->normInf(norm_locs);
    Scalar my_norm_loc = *std::max_element(norm_locs.begin(), norm_locs.end());
    norm = std::max(my_norm_loc, norm);
  }
  return norm;
};


template <typename Scalar>
Scalar
BlockVector<Scalar>::norm1() const
{
  Scalar norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->norm1(norm_locs);
    norm += std::accumulate(norm_locs.begin(), norm_locs.end(), (Scalar)0);
  }
  return norm;
};


template <typename Scalar>
Scalar
BlockVector<Scalar>::norm2() const
{
  Scalar norm = 0.0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> norm_locs(GetComponent_(name)->getNumVectors());
    GetComponent_(name)->norm2(norm_locs());
    norm += std::accumulate(norm_locs.begin(),
                            norm_locs.end(),
                            (Scalar)0,
                            [](Scalar x, Scalar y) { return x + y * y; });
  }
  return sqrt(norm);
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
//     for (int lcv_vector = 0; lcv_vector !=
//     GetComponent_(name)->getNumVectors(); ++lcv_vector) {
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
//     for (int lcv_vector = 0; lcv_vector !=
//     GetComponent_(name)->getNumVectors(); ++lcv_vector) {
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
//     n_loc = GetComponent_(name)->getGlobalLength();
//     for (int lcv_vector = 0; lcv_vector !=
//     GetComponent_(name)->getNumVectors(); ++lcv_vector) {
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
template <typename Scalar>
void
BlockVector<Scalar>::print(std::ostream& os, bool ghosted, bool data_io) const
{
  os << "Comp Vector" << std::endl;
  os << "  components: ";
  for (const auto& name : *this) {
    os << name << "(" << getNumVectors(name) << ") ";
  }
  os << std::endl;

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(os));
  if (data_io) {
    for (const auto& name : *this) { GetComponent_(name, ghosted)->describe(out, Teuchos::VERB_EXTREME); }
  }
};


// Populate by random numbers between -1 and 1.
template <typename Scalar>
void
BlockVector<Scalar>::randomize()
{
  for (const auto& name : *this) { GetComponent_(name)->randomize(); }
};


} // namespace Amanzi

#endif
