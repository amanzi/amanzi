/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  Factory for a CompositeSpace, which is then used to create CompositeVectors.
*/

#ifndef AMANZI_COMPOSITE_VECTOR_SPACE_HH_
#define AMANZI_COMPOSITE_VECTOR_SPACE_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "Mesh.hh"

namespace Amanzi {

class CompositeSpace;
template <typename Scalar>
class CompositeVector_;
using CompositeVector = CompositeVector_<double>;


class CompositeVectorSpace {
 public:
  // Constructor
  CompositeVectorSpace();
  explicit CompositeVectorSpace(const CompositeSpace& map);

  CompositeVectorSpace(const CompositeVectorSpace& other) = default;
  CompositeVectorSpace& operator=(const CompositeVectorSpace&) = default;
  ~CompositeVectorSpace() = default;

  // Checks equality
  bool isSameAs(const CompositeVectorSpace& other) const;
  bool locallySameAs(const CompositeVectorSpace& other) const;
  bool isSubsetOf(const CompositeVectorSpace& other) const;

  //
  // CompositeVectorSpace is primarily a factory.
  // -------------------------------------------------
  // Create a CompositeVector
  template <typename Scalar = double>
  Teuchos::RCP<CompositeVector_<Scalar>> Create() const;

  // Create just the CompositeSpace
  Teuchos::RCP<const CompositeSpace> CreateSpace() const;

  //
  // Specs for the construction of the space
  // -------------------------------------------------
  // Meshes exist on a single communicator.
  Comm_ptr_type getComm() const;

  // mesh specification
  Teuchos::RCP<const AmanziMesh::Mesh> getMesh() const { return mesh_; }
  CompositeVectorSpace* SetMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // ghost spec
  bool Ghosted() const { return ghosted_; }
  CompositeVectorSpace* SetGhosted(bool ghosted = true);

  // Is this space and the resulting CV owned by a PK?
  bool Owned() const { return owned_; }
  CompositeVectorSpace* SetOwned(bool owned = true) {
    owned_ = owned;
    return this;
  }

  // Components are refered to by names.
  // -- Iteration over names of the space
  using name_iterator = std::vector<std::string>::const_iterator;
  name_iterator begin() const { return names_.begin(); }
  name_iterator end() const { return names_.end(); }
  std::size_t size() const { return names_.size(); }

  bool hasComponent(const std::string& name) const
  {
    return indexmap_.find(name) != indexmap_.end();
  }

  // Each component has a number of vectors
  int getNumVectors(const std::string& name) const { return num_vectors_[Index_(name)]; }

  // Each component exists on a mesh entity kind. (CELL, FACE, NODE, etc)
  AmanziMesh::Entity_kind getLocation(const std::string& name) const
  {
    return locations_[Index_(name)];
  }

  BlockMap_ptr_type getMap(const std::string& name, bool ghost = false) const;

  // Update all specs from another space's specs.
  CompositeVectorSpace* Update(const CompositeVectorSpace& other);
  CompositeVectorSpace* Update(const CompositeSpace& other);

  // Update only the components from other
  CompositeVectorSpace* UpdateComponents(const CompositeVectorSpace& other);

  // component specification

  // Add methods append their specs to the space's spec, checking to make
  // sure the spec is OK if the full spec has been set (by an owning PK).
  CompositeVectorSpace*
  AddComponent(const std::string& name, AmanziMesh::Entity_kind location, int num_dofs);

  CompositeVectorSpace* AddComponents(const std::vector<std::string>& names,
                                      const std::vector<AmanziMesh::Entity_kind>& locations,
                                      const std::vector<int>& num_dofs);

  CompositeVectorSpace* AddComponent(const std::string& name,
                                     AmanziMesh::Entity_kind location,
                                     const BlockMap_ptr_type& mastermap,
                                     const BlockMap_ptr_type& ghostmap,
                                     int num_dofs);

  CompositeVectorSpace* AddComponents(const std::vector<std::string>& names,
                                      const std::vector<AmanziMesh::Entity_kind>& locations,
                                      const std::map<std::string, BlockMap_ptr_type>& mastermaps,
                                      const std::map<std::string, BlockMap_ptr_type>& ghostmaps,
                                      const std::vector<int>& num_dofs);

  // Set methods fix the component specs, checking to make sure all previously
  // added specs are contained in the new spec.
  CompositeVectorSpace*
  SetComponent(const std::string& name, AmanziMesh::Entity_kind location, int num_dofs);

  CompositeVectorSpace* SetComponents(const std::vector<std::string>& names,
                                      const std::vector<AmanziMesh::Entity_kind>& locations,
                                      const std::vector<int>& num_dofs);

  CompositeVectorSpace* SetComponent(const std::string& name,
                                     AmanziMesh::Entity_kind location,
                                     const BlockMap_ptr_type& mastermap,
                                     const BlockMap_ptr_type& ghostmap,
                                     int num_dof);

  CompositeVectorSpace* SetComponents(const std::vector<std::string>& names,
                                      const std::vector<AmanziMesh::Entity_kind>& locations,
                                      const std::map<std::string, BlockMap_ptr_type>& mastermaps,
                                      const std::map<std::string, BlockMap_ptr_type>& ghostmaps,
                                      const std::vector<int>& num_dofs);

 private:
  // Indexing of name->std::size_t
  std::size_t Index_(const std::string& name) const
  {
    std::map<std::string, std::size_t>::const_iterator item = indexmap_.find(name);
    AMANZI_ASSERT(item != indexmap_.end());
    return item->second;
  }

  void InitIndexMap_();

  // Consistency checking.
  bool CheckContained_(const std::vector<std::string>& containing,
                       const std::vector<std::string>& contained);

  bool CheckConsistent_(const std::vector<std::string>& names1,
                        const std::vector<AmanziMesh::Entity_kind>& locations1,
                        const std::vector<int>& num_dofs1,
                        const std::vector<std::string>& names2,
                        const std::vector<AmanziMesh::Entity_kind>& locations2,
                        const std::vector<int>& num_dofs2);

  bool UnionAndConsistent_(const std::vector<std::string>& names1,
                           const std::vector<AmanziMesh::Entity_kind>& locations1,
                           const std::vector<int>& num_dofs1,
                           std::vector<std::string>& names2,
                           std::vector<AmanziMesh::Entity_kind>& locations2,
                           std::vector<int>& num_dofs2);

  bool UnionAndConsistent_(const std::vector<std::string>& names1,
                           const std::vector<AmanziMesh::Entity_kind>& locations1,
                           const std::vector<int>& num_dofs1,
                           const std::map<std::string, BlockMap_ptr_type>& mastermaps1,
                           const std::map<std::string, BlockMap_ptr_type>& ghostmaps1,
                           std::vector<std::string>& names2,
                           std::vector<AmanziMesh::Entity_kind>& locations2,
                           std::vector<int>& num_dofs2,
                           std::map<std::string, BlockMap_ptr_type>& mastermaps2,
                           std::map<std::string, BlockMap_ptr_type>& ghostmaps2);

 private:
  bool ghosted_;
  bool owned_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  std::vector<std::string> names_;
  std::map<std::string, std::size_t> indexmap_;

  std::vector<AmanziMesh::Entity_kind> locations_;
  std::vector<int> num_vectors_;

  std::map<std::string, BlockMap_ptr_type> mastermaps_;
  std::map<std::string, BlockMap_ptr_type> ghostmaps_;
};


// Create things.
template <typename Scalar>
Teuchos::RCP<CompositeVector_<Scalar>>
CompositeVectorSpace::Create() const
{
  return Teuchos::rcp(new CompositeVector_<Scalar>(CreateSpace()));
}


} // namespace Amanzi

#endif
