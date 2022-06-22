/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! A patch is data on a collection of meshed entities, as defined by a mesh and a region.

#pragma once

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

#include "Mesh.hh"
#include "CompositeVector.hh"

namespace Amanzi {



//
// A set of entity IDs, defined by a region, mesh, and entity kind.
//
// Note this also has an arbitrary flag_type, which can be a boundary condition
// type or other helper flag.  It is NOT used by this class, but is kept here
// to make using these objects with State easier.
//
struct PatchSpace {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  bool ghosted;
  std::string region;
  AmanziMesh::Entity_kind entity_kind;
  int n_dofs;
  int flag_type;

  PatchSpace() :
    ghosted(false) {}
  PatchSpace(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_,
             bool ghosted_,
             const std::string& region_,
             const AmanziMesh::Entity_kind& entity_kind_,
             const int& n_dofs_,
             const int& flag_type_)
      : mesh(mesh_),
        ghosted(ghosted_),
        region(region_),
        entity_kind(entity_kind_),
        n_dofs(n_dofs_),
        flag_type(flag_type_) {}

  int size() const {
    if (entity_kind == AmanziMesh::BOUNDARY_FACE) {
      Errors::Message msg("Patch cannot handle BOUNDARY_FACE entities, because Mesh does not support sets on these types of entities.  Instead use FACE and filter as needed.");
    }

    return mesh->get_set_size(region, entity_kind,
            ghosted ? AmanziMesh::Parallel_type::ALL : AmanziMesh::Parallel_type::OWNED);
  }
};


//
// A collection of independent patch spaces.  Conceptually these should share
// the same entity_kind, ghosted, mesh, and n_dofs.
struct MultiPatch;

struct MultiPatchSpace {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  bool ghosted;
  int flag_type;
  AmanziMesh::Entity_kind flag_entity;

  MultiPatchSpace() : ghosted(false) {}
  MultiPatchSpace(bool ghosted_) : ghosted(ghosted_) {}
  MultiPatchSpace(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_,
                  bool ghosted_,
                  int flag_type_=-1)
      : mesh(mesh_),
        ghosted(ghosted_),
        flag_type(flag_type_),
        flag_entity(AmanziMesh::UNKNOWN)
  {}

  Teuchos::RCP<MultiPatch> Create() const;

  const PatchSpace& operator[](const int& i) const {
    return subspaces_[i];
  }

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_) {
    mesh = mesh_;
    for (auto& p : *this) {
      p.mesh = mesh_;
    }
  }

  using const_iterator = std::vector<PatchSpace>::const_iterator;
  const_iterator begin() const { return subspaces_.begin(); }
  const_iterator end() const { return subspaces_.end(); }
  std::size_t size() const { return subspaces_.size(); }

  using iterator = std::vector<PatchSpace>::iterator;
  iterator begin() { return subspaces_.begin(); }
  iterator end() { return subspaces_.end(); }

  void AddPatch(const std::string& region,
                AmanziMesh::Entity_kind entity_kind,
                int n_dofs) {
    subspaces_.emplace_back(PatchSpace{mesh, ghosted, region,
                                       entity_kind, n_dofs, flag_type});
  }

 private:
  std::vector<PatchSpace> subspaces_;
};


//
// A set of entity IDs and data on those entities.
//
struct Patch {

  using ViewType = Kokkos::View<double**, Kokkos::LayoutLeft>; 

  Patch(const PatchSpace& space_) :
      space(space_) {
    Kokkos::resize(data, space.size(), space.n_dofs);
  }

  KOKKOS_INLINE_FUNCTION ~Patch(){};

  PatchSpace space;
  // note, this layout is required to ensure that function is slowest-varying,
  // and so can be used with MultiFunction::apply(). See note in
  // MultiFunction.hh
  ViewType data;

  std::size_t size() const { return data.extent(0); }
  std::size_t n_dofs() const { return data.extent(1); }
};

//
// A collection of Patches that share contiguous memory.
//
struct MultiPatch {
  explicit MultiPatch(const MultiPatchSpace& space_) :
      space(space_) {
    for (const auto& subspace : space)
      patches_.emplace_back(Patch{subspace});
  }

  using iterator = typename std::vector<Patch>::iterator;
  iterator begin() { return patches_.begin(); }
  iterator end() { return patches_.end(); }
  std::size_t size() const { return patches_.size(); }

  using const_iterator = typename std::vector<Patch>::const_iterator;
  const_iterator begin() const { return patches_.begin(); }
  const_iterator end() const { return patches_.end(); }

  Patch& operator[](const int& i) {
    return patches_[i];
  }

  const Patch& operator[](const int& i) const {
    return patches_[i];
  }

  MultiPatchSpace space;

 protected:
  std::vector<Patch> patches_;

};




inline Teuchos::RCP<MultiPatch>
MultiPatchSpace::Create() const {
  return Teuchos::rcp(new MultiPatch(*this));
}


inline void
patchToCompositeVector(const Patch& p,
                       const std::string& component,
                       CompositeVector& cv)
{
  auto cv_c = cv.ViewComponent(component, p.space.ghosted);

  // note this is a workaround until we have a way of getting the view on device
  const auto& mesh = cv.Mesh();
  AmanziMesh::Entity_ID_List ids_list;
  mesh->get_set_entities(p.space.region, p.space.entity_kind,
                         p.space.ghosted ? AmanziMesh::Parallel_type::ALL : AmanziMesh::Parallel_type::OWNED,
                         ids_list);

  Kokkos::View<int*,DeviceOnlyMemorySpace> ids("ids", ids_list.size());
  {
    Kokkos::View<int*, MirrorHost, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ids_host(ids_list.data(), ids_list.size());
    Kokkos::deep_copy(ids, ids_host);
  }

  if (component != "boundary_face") {
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>  range({0,0}, {p.data.extent(0), p.data.extent(1)});
    Kokkos::parallel_for(
        "patchToCompositeVector",
        range,
        KOKKOS_LAMBDA(const int& i, const int& j) {
          cv_c(ids(i),j) = p.data(i,j);
        });
  } else {
    AMANZI_ASSERT(false && "Not yet implemented: patchToCompositeVector with boundary_face");
    // have to do some dancing here... this is not correct because p.data is
    // based on faces, but component is based on boundary faces.  Need to
    // either create temporary space, then import, or more likely, unpack the
    // mapping to make sure we only access cv_c on boundary faces.
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>  range({0,0}, {p.data.extent(0), p.data.extent(1)});
    Kokkos::parallel_for(
        "patchToCompositeVector boundary_face",
        range,
        KOKKOS_LAMBDA(const int& i, const int& j) {
          cv_c(ids(i),j) = p.data(i,j);
        });
  }
}

inline void
patchToCompositeVector(const Patch& p,
                       const std::string& component,
                       CompositeVector& cv,
                       CompositeVector_<int>& flag_cv)
{
  // note this is a workaround until we have a way of getting the view on device
  AmanziMesh::Entity_ID_List ids_list;
  cv.Mesh()->get_set_entities(p.space.region, p.space.entity_kind,
                         p.space.ghosted ? AmanziMesh::Parallel_type::ALL : AmanziMesh::Parallel_type::OWNED,
                         ids_list);

  Kokkos::View<int*,DeviceOnlyMemorySpace> ids("ids", ids_list.size());
  {
    Kokkos::View<int*, MirrorHost, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ids_host(ids_list.data(), ids_list.size());
    Kokkos::deep_copy(ids, ids_host);
  }

  if (component != "boundary_face") {
    AMANZI_ASSERT(ids.extent(0) == p.data.extent(0));
    auto flag_type = p.space.flag_type;

    auto cv_c = cv.ViewComponent(component, p.space.ghosted);
    auto flag_c = flag_cv.ViewComponent(component, p.space.ghosted);

    Kokkos::parallel_for(
        "patchToCompositeVector",
        ids.extent(0),
        KOKKOS_LAMBDA(const int& i) {
          cv_c(ids(i),0) = p.data(i,0);
          flag_c(ids(i),0) = flag_type;
        });
  } else {
    AMANZI_ASSERT(false && "Not yet implemented: patchToCompositeVector with boundary_face");
    // have to do some dancing here... this is not correct because p.data is
    // based on faces, but component is based on boundary faces.  Need to
    // either create temporary space, then import, or more likely, unpack the
    // mapping to make sure we only access cv_c on boundary faces.
    // Kokkos::parallel_for(
    //     "patchToCompositeVector boundary_face",
    //     p.data.extent(0),
    //     KOKKOS_LAMBDA(const int& i) {
    //       cv_c(ids(i),0) = p.data(i,0);
    //       flag_c(ids(i), 0) = p.space.flag_type;
    //     });
  }
}


//
// Copies values from a set of patches into a vector.
//
inline void
multiPatchToCompositeVector(const MultiPatch& mp,
                            const std::string& component,
                            CompositeVector& cv)
{
  for (const auto& p : mp) {
    patchToCompositeVector(p, component, cv);
  }
}

//
// Copies values and flag from a set of patches into a vector and a flag vector.
//
void inline
multiPatchToCompositeVector(const MultiPatch& mp,
                            const std::string& component,
                            CompositeVector& cv,
                            CompositeVector_<int>& flag)
{
  for (const auto& p : mp) {
    patchToCompositeVector(p, component, cv, flag);
  }
}

} // namespace Amanzi
