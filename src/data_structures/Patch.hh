/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! A patch is data on a collection of meshed entities, as defined by a mesh and
//! a region.
#pragma once

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

#include "Key.hh"
#include "Mesh.hh"

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
  int num_vectors;
  int flag_type;
  Kokkos::View<int*> flags;


  PatchSpace() : ghosted(false) {}
  PatchSpace(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_,
             bool ghosted_,
             const std::string& region_,
             const AmanziMesh::Entity_kind& entity_kind_,
             const int& num_vectors_,
             const int& flag_type_)
    : mesh(mesh_),
      ghosted(ghosted_),
      region(region_),
      entity_kind(entity_kind_),
      num_vectors(num_vectors_),
      flag_type(flag_type_)
  {
    if (flag_type == -1) Kokkos::resize(flags, size());
  }
  PatchSpace(const PatchSpace& other) = default;
  PatchSpace& operator=(const PatchSpace&) = default;
  ~PatchSpace() = default;

  bool operator==(const PatchSpace& other) const
  {
    return mesh == other.mesh && region == other.region && entity_kind == other.entity_kind &&
           num_vectors == other.num_vectors;
  }

  int size() const
  {
    if (entity_kind == AmanziMesh::BOUNDARY_FACE) {
      Errors::Message msg("Patch cannot handle BOUNDARY_FACE entities, because "
                          "Mesh does not support sets on these types of "
                          "entities.  Instead use FACE and filter as needed.");
    }

    return mesh->getSetSize(region,
                            entity_kind,
                            ghosted ? AmanziMesh::Parallel_kind::ALL :
                                      AmanziMesh::Parallel_kind::OWNED);
  }

  AmanziMesh::Mesh::cEntity_ID_View getIDs() const
  {
    return mesh->getSetEntities(region,
                                entity_kind,
                                ghosted ? AmanziMesh::Parallel_kind::ALL :
                                          AmanziMesh::Parallel_kind::OWNED);
  }
};


//
// A collection of independent patch spaces.  Conceptually these may share the
// same entity_kind, ghosted, mesh, and num_vectors, or may not, but they MUST
// share the same mesh.
//
template <typename T>
struct MultiPatch;

struct MultiPatchSpace {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  AmanziMesh::Entity_kind entity_kind;
  bool ghosted;
  int flag_type;

  MultiPatchSpace() : ghosted(false), entity_kind(AmanziMesh::Entity_kind::UNKNOWN) {}
  MultiPatchSpace(bool ghosted_) : ghosted(ghosted_), entity_kind(AmanziMesh::Entity_kind::UNKNOWN)
  {}
  MultiPatchSpace(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_,
                  bool ghosted_,
                  AmanziMesh::Entity_kind entity_kind_ = AmanziMesh::Entity_kind::UNKNOWN,
                  int flag_type_ = -1)
    : mesh(mesh_), entity_kind(entity_kind_), ghosted(ghosted_), flag_type(flag_type_)
  {}

  Teuchos::RCP<const PatchSpace> operator[](const int& i) const { return subspaces_[i]; }

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_)
  {
    mesh = mesh_;
    for (auto& p : *this) { p->mesh = mesh_; }
  }
  void set_flag(int flag)
  {
    flag_type = flag;
    for (auto& p : *this) { p->flag_type = flag; }
  }
  void set_entity_kind(AmanziMesh::Entity_kind entity_kind_)
  {
    entity_kind = entity_kind_;
    for (auto& p : *this) { p->entity_kind = entity_kind_; }
  }

  using const_iterator = std::vector<Teuchos::RCP<PatchSpace>>::const_iterator;
  const_iterator begin() const { return subspaces_.begin(); }
  const_iterator end() const { return subspaces_.end(); }
  std::size_t size() const { return subspaces_.size(); }

  using iterator = std::vector<Teuchos::RCP<PatchSpace>>::const_iterator;
  iterator begin() { return subspaces_.begin(); }
  iterator end() { return subspaces_.end(); }

  void addPatch(const std::string& region, AmanziMesh::Entity_kind entity_kind, int num_vectors)
  {
    subspaces_.emplace_back(
      Teuchos::rcp(new PatchSpace{ mesh, ghosted, region, entity_kind, num_vectors, flag_type }));
  }

  void addPatch(const Teuchos::RCP<PatchSpace>& ps) { subspaces_.emplace_back(ps); }

 private:
  std::vector<Teuchos::RCP<PatchSpace>> subspaces_;
};


//
// A set of entity IDs and data on those entities.
//
template <typename T>
struct Patch {
  using View_type = Kokkos::View<T**, Kokkos::LayoutLeft>;

  Patch(const Teuchos::RCP<const PatchSpace>& space_) : space(space_)
  {
    Kokkos::resize(data, space->size(), space->num_vectors);
  }

  Patch(const Patch& other) = default;
  Patch& operator=(const Patch&) = default;
  ~Patch() = default;

  Teuchos::RCP<const PatchSpace> space;
  // note, this layout is required to ensure that function is slowest-varying,
  // and so can be used with MultiFunction::apply(). See note in
  // MultiFunction.hh
  View_type data;

  std::size_t size() const { return data.extent(0); }
  std::size_t getNumVectors() const { return data.extent(1); }
};

//
// A collection of Patches that share contiguous memory.
//
template <typename T>
struct MultiPatch {
  explicit MultiPatch(const Teuchos::RCP<const MultiPatchSpace>& space_) : space(space_)
  {
    for (const auto& subspace : *space) patches_.emplace_back(Patch<T>{ subspace });
  }

  using iterator = typename std::vector<Patch<T>>::iterator;
  iterator begin() { return patches_.begin(); }
  iterator end() { return patches_.end(); }
  std::size_t size() const { return patches_.size(); }

  using const_iterator = typename std::vector<Patch<T>>::const_iterator;
  const_iterator begin() const { return patches_.begin(); }
  const_iterator end() const { return patches_.end(); }

  Patch<T>& operator[](const int& i) { return patches_[i]; }

  const Patch<T>& operator[](const int& i) const { return patches_[i]; }

  Teuchos::RCP<const MultiPatchSpace> space;

 protected:
  std::vector<Patch<T>> patches_;
};


//
// Nonmember function for creating PatchSpaces from a ParameterList
//
inline void
readMultiPatchSpace(Teuchos::ParameterList& list, MultiPatchSpace& mps)
{
  Teuchos::Array<std::string> regions;
  if (list.isParameter("regions")) {
    regions = list.get<Teuchos::Array<std::string>>("regions");
  } else {
    regions.push_back(list.get<std::string>("region", Keys::cleanPListName(list.name())));
  }

  Teuchos::Array<AmanziMesh::Entity_kind> entity_kinds;
  if (mps.entity_kind == AmanziMesh::Entity_kind::UNKNOWN) {
    if (list.isParameter("entity kind")) {
      entity_kinds.push_back(AmanziMesh::createEntityKind(list.get<std::string>("entity kind")));
    } else {
      auto ekinds = list.get<Teuchos::Array<std::string>>("entity kinds");
      for (auto ekind : ekinds) entity_kinds.push_back(AmanziMesh::createEntityKind(ekind));
    }
  } else {
    entity_kinds.push_back(mps.entity_kind);
  }

  int n_dofs = list.get<int>("number of vectors", 1);

  for (auto r : regions) {
    for (auto e_kind : entity_kinds) { mps.addPatch(r, e_kind, n_dofs); }
  }
}

//
// Nonmember function for creating MultiPatchSpaces from a ParameterList
//
inline Teuchos::RCP<MultiPatchSpace>
createMultiPatchSpace(Teuchos::ParameterList& list,
                      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                      AmanziMesh::Entity_kind entity_kind = AmanziMesh::Entity_kind::UNKNOWN,
                      int flag_type = 0)
{
  // All are expected to be sublists of identical structure.
  auto mps = Teuchos::rcp(new MultiPatchSpace(mesh, false, entity_kind, flag_type));

  for (auto sublist : list) {
    std::string name = sublist.first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec_plist = list.sublist(name);

      try {
        readMultiPatchSpace(spec_plist, *mps);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist " << name << ": " << msg.what();
        throw(m);
      }

    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter " << name << " is not a sublist";
      throw(m);
    }
  }
  return mps;
}


} // namespace Amanzi
