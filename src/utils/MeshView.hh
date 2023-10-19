/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once

#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

#include "Iterators.hh"

namespace Kokkos {

template <class DataType, class... Properties>
struct MeshView : public Kokkos::View<DataType, Properties...> {
  using baseView = Kokkos::View<DataType, Properties...>;
  using baseView::baseView;
  using iterator = Impl::View_iter<MeshView<DataType, Properties...>>;
  using traits = ViewTraits<DataType, Properties...>;

  using const_type = MeshView<typename traits::const_data_type,
                              typename traits::array_layout,
                              typename traits::device_type,
                              typename traits::memory_traits>;
  using const_iterator = Impl::View_iter<const_type>;

  MeshView(const baseView& bv) : baseView(bv) {}
  MeshView(const MeshView& bv) : baseView(bv) {}

  KOKKOS_FUNCTION MeshView& operator=(const MeshView& other)
  {
    baseView::operator=(other);
    return *this;
  }

  KOKKOS_FUNCTION MeshView& operator=(const MeshView&& other)
  {
    baseView::operator=(other);
    return *this;
  }

  using HostMirror =
    MeshView<typename traits::non_const_data_type,
             typename traits::array_layout,
             Device<DefaultHostExecutionSpace, typename traits::host_mirror_space::memory_space>>;

  KOKKOS_INLINE_FUNCTION iterator begin() const { return iterator(*this); }

  KOKKOS_INLINE_FUNCTION iterator end() const { return iterator(*this, this->size()); }

  KOKKOS_INLINE_FUNCTION const_iterator cbegin() const { return const_iterator(*this); }

  KOKKOS_INLINE_FUNCTION const_iterator cend() const { return const_iterator(*this, this->size()); }

  void insert(iterator v0_e, iterator v1_b, iterator v1_e)
  {
    //assert(v0_e - *this->end() != 0 && "Only insert at end supported for MeshViews");
    std::size_t size = v1_e - v1_b;
    std::size_t csize = this->size();
    Kokkos::resize(*this, size + csize);
    for (int i = csize; i < this->size(); ++i, ++v1_b) this->operator[](i) = *(v1_b);
  }

  template <typename MV>
  KOKKOS_INLINE_FUNCTION void fromConst(const MV& cmv)
  {
    Kokkos::resize(*this, cmv.size());
    Kokkos::deep_copy(*this, cmv);
  }
};

template <class SubDataType, class... SubProperties, class... Args>
MeshView<SubDataType, SubProperties...>
subview(Kokkos::MeshView<SubDataType, SubProperties...> v, Args... args)
{
  MeshView ret(
    Kokkos::subview((typename Kokkos::MeshView<SubDataType, SubProperties...>::baseView)v,
                    std::forward<Args>(args)...));
  return ret;
}

template <class DataType, class Arg1Type = void, class Arg2Type = void, class Arg3Type = void>
struct MeshDualView : public Kokkos::ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type> {
  using traits = Kokkos::ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type>;
  using host_mirror_space = typename traits::host_mirror_space;
  using t_dev = MeshView<typename traits::data_type, Arg2Type, Arg2Type, Arg3Type>;
  using t_host = typename t_dev::HostMirror;
  using t_modified_flags = MeshView<unsigned int[2], LayoutLeft, Kokkos::HostSpace>;
  t_modified_flags modified_flags;

  KOKKOS_INLINE_FUNCTION t_host view_host() const { return h_view; }
  KOKKOS_INLINE_FUNCTION t_dev view_device() const { return d_view; }

  t_dev d_view;
  t_host h_view;

  MeshDualView() = default;
  MeshDualView(const std::string& label, const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
    : modified_flags(t_modified_flags("DualView::modified_flags")),
      d_view(label, n0),
      h_view(create_mirror_view(d_view)) // without UVM, host View mirrors
  {}
  MeshDualView(const MeshDualView& src)
    : modified_flags(src.modified_flags), d_view(src.d_view), h_view(src.h_view)
  {}
  template <class SS, class LS, class DS, class MS>
  MeshDualView(const MeshDualView<SS, LS, DS, MS>& src)
    : modified_flags(src.modified_flags), d_view(src.d_view), h_view(src.h_view)
  {}

  void resize(const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
  {
    if (modified_flags.data() == nullptr) {
      modified_flags = t_modified_flags("DualView::modified_flags");
    }
    if (modified_flags(1) >= modified_flags(0)) {
      /* Resize on Device */
      ::Kokkos::resize(d_view, n0);
      h_view = create_mirror_view(d_view);

      /* Mark Device copy as modified */
      modified_flags(1) = modified_flags(1) + 1;

    } else {
      /* Realloc on Device */

      ::Kokkos::realloc(d_view, n0);

      const bool sizeMismatch = (h_view.extent(0) != n0);
      if (sizeMismatch) ::Kokkos::resize(h_view, n0);

      t_host temp_view = create_mirror_view(d_view);

      /* Remap on Host */
      Kokkos::deep_copy(temp_view, h_view);

      h_view = temp_view;

      d_view = create_mirror_view(typename t_dev::execution_space(), h_view);

      /* Mark Host copy as modified */
      modified_flags(0) = modified_flags(0) + 1;
    }
  }
};

} // namespace Kokkos
