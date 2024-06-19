/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
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

#include "AmanziMap.hh"
#include "AmanziVector.hh"


namespace Amanzi {

// Constructor
template <typename Scalar>
BlockVector<Scalar>::BlockVector(const Teuchos::RCP<const BlockSpace>& map, InitMode mode)
  : map_(map)
{
  for (const auto& name : *this) {
    setComponent_(name, true, Teuchos::null);
    setComponent_(name, false, Teuchos::null);
  }
  CreateData_(mode);
};


// copy constructor
template <typename Scalar>
BlockVector<Scalar>::BlockVector(const BlockVector<Scalar>& other,
                                 Teuchos::DataAccess access,
                                 InitMode mode)
  : map_(other.getMap())
{
  for (const auto& name : *this) {
    setComponent_(name, true, Teuchos::null);
    setComponent_(name, false, Teuchos::null);
  }

  if (access == Teuchos::DataAccess::View) {
    Errors::Message message("BlockVector: View semantic not supported.");
    throw(message);

    // for (const auto& name : *this) {
    //   setComponent_(name, true, other.getComponent_(name, true));
    //   setComponent_(name, false, other.getComponent_(name, false));
    // }
  } else {
    CreateData_(mode);
    if (mode == InitMode::COPY) { *this = other; }
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
      auto ghost_map = getMap()->getComponentMap(name, true);
      if (ghost_map != Teuchos::null) {
        setComponent_(
          name,
          true,
          Teuchos::rcp(new MultiVector_type_<Scalar>(ghost_map, getNumVectors(name), true)));

        auto g_comp = getComponent_(name, true);
        auto master_map = getMap()->getComponentMap(name, false);

        // get the ghost component's data
        auto m_comp = g_comp->offsetViewNonConst(master_map, 0);
        setComponent_(name, false, m_comp);

      } else {
        auto master_map = getMap()->getComponentMap(name, false);
        setComponent_(
          name,
          false,
          Teuchos::rcp(new MultiVector_type_<Scalar>(master_map, getNumVectors(name), true)));
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
      const auto other_v = other.getComponent(name, false);
      const auto this_v = getComponent(name, false);
      if (other_v == Teuchos::null && this_v == Teuchos::null) {
        // pass
      } else if (this_v == Teuchos::null) {
        auto new_v = Teuchos::rcp(new MultiVector_type_<Scalar>(*other_v, Teuchos::Copy));
        Tpetra::deep_copy(*new_v, *other_v);
        setComponent_(name, false, new_v);
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
Teuchos::RCP<const MultiVector_type_<Scalar>>
BlockVector<Scalar>::getComponent(const std::string& name, bool ghosted) const
{
  if (!getMap()->hasComponent(name)) {
    Errors::Message message("BlockVector: Requested component (" + name + ") does not exist.");
    throw(message);
  }
  return getComponent_(name, ghosted);
}


// View data, non-const version.
template <typename Scalar>
Teuchos::RCP<MultiVector_type_<Scalar>>
BlockVector<Scalar>::getComponent(const std::string& name, bool ghosted)
{
  if (!getMap()->hasComponent(name)) {
    Errors::Message message("BlockVector: Requested component (" + name + ") does not exist.");
    throw(message);
  }
  return getComponent_(name, ghosted);
};


// // Set data
// template<typename Scalar>
// void
// BlockVector<Scalar>::setComponent(const std::string& name,
//         const Teuchos::RCP<MultiVector_type_<Scalar>>& data)
// {
//   if (!getMap()->hasComponent(name)) {
//     Errors::Message message("BlockVector: Requested setComponent("+name+")
//     does not exist."); throw(message);
//   }
//   if (data != Teuchos::null) {
//     if (!getComponentMap(name)->isSameAs(*data->getMap())) {
//       Errors::Message message("BlockVector: Requested setComponent("+name+")
//       has incompatible map."); throw(message);
//     }
//     if (!(getNumVectors(name) == data->getNumVectors())) {
//       Errors::Message message("BlockVector: Requested setComponent("+name+")
//       has incompatible number of DoFs."); throw(message);
//     }
//   }
//   data_[name] = data;
// }

// -- View a component vector
template <typename Scalar>
template<MemSpace_kind MEM>
auto
BlockVector<Scalar>::viewComponent(const std::string& name, bool ghosted) const
{
  if constexpr(MEM == MemSpace_kind::HOST) {
    return getComponent_(name, ghosted)
      ->template getLocalView<typename HostView_type::memory_space>(Tpetra::Access::ReadOnly);
  } else {
    return getComponent_(name, ghosted)
      ->template getLocalView<typename View_type::memory_space>(Tpetra::Access::ReadOnly);
  }
}

template <typename Scalar>
template<MemSpace_kind MEM>
auto
BlockVector<Scalar>::viewComponent(const std::string& name, bool ghosted)
{
  if constexpr(MEM == MemSpace_kind::HOST) {
    return getComponent_(name, ghosted)
      ->template getLocalView<typename HostView_type::memory_space>(Tpetra::Access::ReadWrite);
  } else {
    return getComponent_(name, ghosted)
      ->template getLocalView<typename View_type::memory_space>(Tpetra::Access::ReadWrite);
  }
}

// -- SubView of a component vector
template <typename Scalar>
template<MemSpace_kind MEM>
auto
BlockVector<Scalar>::viewComponent(const std::string& name, std::size_t dof, bool ghosted) const
{
  return Kokkos::subview(viewComponent<MEM>(name, ghosted), Kokkos::ALL(), dof);
}

template <typename Scalar>
template<MemSpace_kind MEM>
auto
BlockVector<Scalar>::viewComponent(const std::string& name, std::size_t dof, bool ghosted)
{
  return Kokkos::subview(viewComponent<MEM>(name, ghosted), Kokkos::ALL(), dof);
}


// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
template <typename Scalar>
void
BlockVector<Scalar>::scatterMasterToGhosted(Tpetra::CombineMode mode) const
{
  for (const auto& name : *this) { scatterMasterToGhosted(name, mode); }
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
BlockVector<Scalar>::scatterMasterToGhosted(const std::string& name, Tpetra::CombineMode mode) const
{
  auto import = getMap()->getImporter(name);
  if (import == Teuchos::null) return;
  auto g_comp = ghost_data_.at(name); // note cannot use getComponent_() here, need non-const
  if (g_comp == Teuchos::null) return;
  auto m_comp = getComponent_(name, false);
  g_comp->doImport(*m_comp, *import, mode);
};


// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is Add,
// where off-process values are first summed into the on-process value.
template <typename Scalar>
void
BlockVector<Scalar>::gatherGhostedToMaster(Tpetra::CombineMode mode)
{
  for (const auto& name : *this) { gatherGhostedToMaster(name, mode); }
};


template <typename Scalar>
void
BlockVector<Scalar>::gatherGhostedToMaster(const std::string& name, Tpetra::CombineMode mode)
{
  auto import = getMap()->getImporter(name);
  if (import == Teuchos::null) return;
  // communicate
  auto g_comp = getComponent(name, true);
  auto m_comp = getComponent(name, false);
  m_comp->doExport(*g_comp, *import, mode);
};


// Vector operations.
// -- Insert value into data.
template <typename Scalar>
void
BlockVector<Scalar>::putScalar(Scalar scalar)
{
  for (const auto& name : *this) { getComponent_(name)->putScalar(scalar); }
};


// -- Insert value into component [name].
template <typename Scalar>
void
BlockVector<Scalar>::putScalar(const std::string& name, Scalar scalar)
{
  getComponent(name)->putScalar(scalar);
};


// -- Insert values into data of component [name].
template <typename Scalar>
void
BlockVector<Scalar>::putScalar(const std::string& name, const std::vector<Scalar>& scalar)
{
  if (scalar.size() != getNumVectors(name)) {
    Errors::Message message("BlockVector: requested putScalar(" + name +
                            ") with incorrectly sized vector of scalars.");
    throw(message);
  }

  for (int lcv_vector = 0; lcv_vector != getNumVectors(name); ++lcv_vector) {
    (*getComponent_(name))(lcv_vector)->putScalar(scalar[lcv_vector]);
  }
};


// Vector operations.
// -- Insert value into data.
template <typename Scalar>
void
BlockVector<Scalar>::putScalarMasterAndGhosted(Scalar scalar)
{
  for (const auto& name : *this) { getComponent_(name, true)->putScalar(scalar); }
};


template <typename Scalar>
void
BlockVector<Scalar>::putScalarGhosted(Scalar scalar)
{
  for (const auto& comp : *this) {
    auto vv = viewComponent(comp, true);

    int size_owned = getComponent(comp, false)->getLocalLength();
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> range_policy({ size_owned, 0 },
                                                        { vv.extent(0), vv.extent(1) });
    Kokkos::parallel_for(
      "BlockVector_impl::putScalarGhosted", range_policy, KOKKOS_LAMBDA(const int i, const int j) {
        vv(i, j) = scalar;
      });
  }
}


// this <- abs(this)
template <typename Scalar>
void
BlockVector<Scalar>::abs(const BlockVector<Scalar>& other)
{
  for (const auto& name : *this) { getComponent_(name)->abs(*other.getComponent_(name)); }
}


// -- this <- value*this
template <typename Scalar>
void
BlockVector<Scalar>::scale(Scalar value)
{
  for (const auto& name : *this) { getComponent_(name)->scale(value); }
};


// scale() applied to component name.
template <typename Scalar>
void
BlockVector<Scalar>::scale(const std::string& name, Scalar value)
{
  getComponent_(name)->scale(value);
};


// this <- abs(this)
template <typename Scalar>
void
BlockVector<Scalar>::reciprocal(const BlockVector<Scalar>& other)
{
  for (const auto& name : *this) { getComponent_(name)->reciprocal(*other.getComponent_(name)); }
}


// -- result <- other \dot this
template <typename Scalar>
Scalar
BlockVector<Scalar>::dot(const BlockVector<Scalar>& other) const
{
  Scalar result = 0;
  for (const auto& name : *this) {
    Teuchos::Array<Scalar> intermediate_result(getNumVectors(name));
    getComponent_(name)->dot(*(other.getComponent_(name)), intermediate_result);
    result += std::accumulate(intermediate_result.begin(), intermediate_result.end(), (Scalar)0);
  }
  return result;
};


// -- this <- scalarA*A + scalarThis*this
template <typename Scalar>
void
BlockVector<Scalar>::update(Scalar scalarA, const BlockVector<Scalar>& A, Scalar scalarThis)
{
  for (const auto& name : *this) {
    if (A.hasComponent(name)) {
      getComponent_(name)->update(scalarA, *A.getComponent_(name), scalarThis);
    }
  }
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
template <typename Scalar>
void
BlockVector<Scalar>::update(Scalar scalarA,
                            const BlockVector<Scalar>& A,
                            Scalar scalarB,
                            const BlockVector<Scalar>& B,
                            Scalar scalarThis)
{
  for (const auto& name : *this) {
    if (A.hasComponent(name) && B.hasComponent(name)) {
      getComponent_(name)->update(
        scalarA, *A.getComponent_(name), scalarB, *B.getComponent_(name), scalarThis);
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
    if (B.getComponent_(name) != this->getComponent(name) &&
        A.getComponent_(name) == this->getComponent(name)) {
      if (B.getComponent_(name)->getNumVectors() > 1) {
        Errors::Message message("Not implemented multiply: Tpetra does not "
                                "provide multivector elementwise multiply.");
        Exceptions::amanzi_throw(message);
      }

      getComponent_(name)->elementWiseMultiply(
        scalarAB, *B.getComponent_(name)->getVector(0), *A.getComponent_(name), scalarThis);
    } else if (B.getComponent_(name) == this->getComponent(name) &&
               A.getComponent_(name) == this->getComponent(name)) {
      Errors::Message message("Unclear whether, in elementWiseMultiply(), "
                              "A,B,& this can all be the same vector.");
      Exceptions::amanzi_throw(message);
    } else {
      if (A.getComponent_(name)->getNumVectors() > 1) {
        Errors::Message message("Not implemented multiply: Tpetra does not "
                                "provide multivector elementwise multiply.");
        Exceptions::amanzi_throw(message);
      }

      getComponent_(name)->elementWiseMultiply(
        scalarAB, *A.getComponent_(name)->getVector(0), *B.getComponent_(name), scalarThis);
    }
  }
};

// // -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise
// division template<typename Scalar> int
// BlockVector<Scalar>::ReciprocalelementWiseMultiply(Scalar scalarAB, const
// BlockVector<Scalar>& A, const BlockVector<Scalar>& B,
//                   Scalar scalarThis) {
//   for (const auto& name : *this) {
//     getComponent_(name)->reciprocalelementWiseMultiply(scalarAB,
//     *A.getComponent_(name), *B.getComponent_(name), scalarThis);
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
    Teuchos::Array<Scalar> norm_locs(getComponent_(name)->getNumVectors());
    getComponent_(name)->normInf(norm_locs);
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
    Teuchos::Array<Scalar> norm_locs(getComponent_(name)->getNumVectors());
    getComponent_(name)->norm1(norm_locs);
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
    Teuchos::Array<Scalar> norm_locs(getComponent_(name)->getNumVectors());
    getComponent_(name)->norm2(norm_locs());
    norm += std::accumulate(
      norm_locs.begin(), norm_locs.end(), (Scalar)0, [](Scalar x, Scalar y) { return x + y * y; });
  }
  return sqrt(norm);
};


// Debugging?
template <typename Scalar>
void
BlockVector<Scalar>::print(std::ostream& os, bool ghosted, bool data_io) const
{
  os << "Comp Vector" << std::endl;
  os << "  components: ";
  for (const auto& name : *this) { os << name << "(" << getNumVectors(name) << ") "; }
  os << std::endl;

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(os));
  if (data_io) {
    for (const auto& name : *this) {
      getComponent_(name, ghosted)->describe(out, Teuchos::VERB_EXTREME);
    }
  }
};


// Populate by random numbers between -1 and 1.
template <typename Scalar>
void
BlockVector<Scalar>::randomize()
{
  for (const auto& name : *this) { getComponent_(name)->randomize(); }
};


} // namespace Amanzi

#endif
