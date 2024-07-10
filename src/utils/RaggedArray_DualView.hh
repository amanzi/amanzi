/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once

#include "MeshView.hh"
#include "ViewUtils.hh"

namespace Amanzi {

template <typename T>
struct RaggedArray_DualView {
  using View_type = Kokkos::MeshView<T*, Kokkos::DefaultHostExecutionSpace>;
  using Entity_ID_View = View_type;

  using type_t = T;

  template <MemSpace_kind MEM>
  using constview = Kokkos::MeshView<const T*, MemoryLocation<MEM>>;

  Kokkos::MeshDualView<int*> rows;
  Kokkos::MeshDualView<T*> entries;

  using host_mirror_space = typename Kokkos::MeshDualView<T*>::host_mirror_space;
  using execution_space = typename Kokkos::MeshDualView<T*>::execution_space;

  RaggedArray_DualView()
    : rows("RaggedArray_DualView: Rows", 0), entries("RaggedArray_DualView: Entries", 0)
  {}


  RaggedArray_DualView(const std::vector<std::vector<T>>& vect)
  {
    rows.resize(vect.size() + 1);
    int count = 0;
    view<MemSpace_kind::HOST>(rows)[0] = count;
    for (int i = 1; i < vect.size() + 1; ++i) {
      count += vect[i - 1].size();
      view<MemSpace_kind::HOST>(rows)[i] = count;
    }

    entries.resize(count);

    int cur = 0;
    for (int i = 0; i < vect.size(); ++i) {
      for (int j = 0; j < vect[i].size(); ++j) {
        view<MemSpace_kind::HOST>(entries)[cur++] = vect[i][j];
      }
    }
    update<MemSpace_kind::DEVICE>();
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRowUnmanaged(int row)
  {
    const std::size_t size = view<MEM>(rows)[row + 1] - view<MEM>(rows)[row];
    const auto ptr = view<MEM>(entries).data() + view<MEM>(rows)[row];
    return Kokkos::MeshView<T*, MemoryLocation<MEM>>(ptr, size);
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRowUnmanaged(int row) const
  {
    const std::size_t size = view<MEM>(rows)[row + 1] - view<MEM>(rows)[row];
    const auto ptr = view<MEM>(entries).data() + view<MEM>(rows)[row];
    return Kokkos::MeshView<T*, MemoryLocation<MEM>>(ptr, size);
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRow(int row)
  {
    return Kokkos::subview(view<MEM>(entries),
                           Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row + 1]));
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRow(int row) const
  {
    return Kokkos::subview(view<MEM>(entries),
                           Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row + 1]));
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION T& get(int row, int i)
  {
    return view<MEM>(entries)[view<MEM>(rows)[row] + i];
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION const T& get(int row, int i) const
  {
    return view<MEM>(entries)[view<MEM>(rows)[row] + i];
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION int size() const
  {
    return view<MEM>(rows).size() - 1;
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION int size(int row) const
  {
    return view<MEM>(rows)[row + 1] - view<MEM>(rows)[row];
  }

  template <MemSpace_kind MEM>
  void update()
  {
    if constexpr (MEM == MemSpace_kind::HOST) {
      Kokkos::deep_copy(rows.view_host(), rows.view_device());
      Kokkos::deep_copy(entries.view_host(), entries.view_device());
    } else {
      Kokkos::deep_copy(rows.view_device(), rows.view_host());
      Kokkos::deep_copy(entries.view_device(), entries.view_host());
    }
  }

  bool operator!=(const RaggedArray_DualView& oth)
  {
    return oth.rows.view_host().size() != rows.view_host().size() &&
           oth.entries.view_host().size() != entries.view_host().size();
  }

  void resize(int s1, int s2)
  {
    rows.resize(s1 + 1);
    entries.resize(s1 * s2);
    rows.h_view[0] = 0;
    for (int i = 1; i < s1 + 1; ++i) { rows.h_view[i] = rows.h_view[i - 1] + s2; }
    update<MemSpace_kind::DEVICE>();
  }
};

//
// Cache a RaggedArray from a callable, e.g. getCellFaces()
//
template <typename T, typename Func>
auto
asRaggedArray_DualView(Func mesh_func, int count)
{
  RaggedArray_DualView<T> adj;
  adj.rows.resize(count + 1);

  // do a count first, setting rows
  std::vector<Kokkos::MeshView<const T*, Kokkos::HostSpace>> ents(count);
  int total = 0;
  for (int i = 0; i != count; ++i) {
    view<MemSpace_kind::HOST>(adj.rows)[i] = total;

    mesh_func(i, ents[i]);
    total += ents[i].size();
  }
  view<MemSpace_kind::HOST>(adj.rows)[count] = total;
  adj.entries.resize(total);

  for (int i = 0; i != count; ++i) {
    const auto& ent = ents[i];
    for (int j = 0; j != ent.size(); ++j) { adj.template get<MemSpace_kind::HOST>(i, j) = ent[j]; }
  }

  Kokkos::deep_copy(adj.rows.view_device(), adj.rows.view_host());
  Kokkos::deep_copy(adj.entries.view_device(), adj.entries.view_host());
  return adj;
}

} // namespace Amanzi
