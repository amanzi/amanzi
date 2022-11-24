/*
  Data Structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Daniil Svyatsky (dasvyat@lanl.gov)

  Factory for a CompositeVector on an Amanzi mesh.

  This should be thought of as a vector-space: it lays out data components as a
  mesh along with entities on the mesh.  This meta-data can be used with the
  mesh's *_map() methods to create the data.

  This class is very light weight as it maintains only meta-data.
*/

#ifndef AMANZI_COMPOSITEVECTOR_SPACE_HH_
#define AMANZI_COMPOSITEVECTOR_SPACE_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "Mesh.hh"

namespace Amanzi {

class CompositeVector;

// Nonmember helper function
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map>>
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location);


class CompositeVectorSpace {
 public:
  // Constructor
  CompositeVectorSpace();
  CompositeVectorSpace(const CompositeVectorSpace& other) = default;
  CompositeVectorSpace(const CompositeVectorSpace& other, bool ghosted);

  // assignment
  CompositeVectorSpace& operator=(const CompositeVectorSpace&) = default;

  // CompositeVectorSpace is a factory of CompositeVectors
  Teuchos::RCP<CompositeVector> Create() const;

  // Checks equality
  bool SameAs(const CompositeVectorSpace& other) const;
  bool SubsetOf(const CompositeVectorSpace& other) const;

  // -------------------------------------
  // Specs for the construction of CVs
  // -------------------------------------

  // CompositeVectors exist on a single communicator.
  Comm_ptr_type Comm() const { return mesh_->get_comm(); }

  // mesh specification
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return mesh_; }
  CompositeVectorSpace* SetMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // ghost spec
  bool Ghosted() const { return ghosted_; }
  CompositeVectorSpace* SetGhosted(bool ghosted = true);

  // Is this space and the resulting CV owned by a PK?
  bool Owned() const { return owned_; }
  void SetOwned(bool owned = true) { owned_ = owned; }

  // Components are refered to by names.
  // -- Iteration over names of the space
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return names_.begin(); }
  name_iterator end() const { return names_.end(); }
  unsigned int size() const { return names_.size(); }

  bool HasComponent(const std::string& name) const
  {
    return indexmap_.find(name) != indexmap_.end();
  }
  int NumComponents() const { return size(); }

  // Each component has a number of Degrees of Freedom.
  int NumVectors(const std::string& name) const { return num_dofs_[Index_(name)]; }

  // Each component exists on a mesh entity kind. (CELL, FACE, NODE, etc)
  AmanziMesh::Entity_kind Location(const std::string& name) const
  {
    return locations_[Index_(name)];
  }

  Teuchos::RCP<const Epetra_BlockMap> Map(const std::string& name, bool ghost = false) const;

  // Update all specs from another space's specs.
  // Useful for PKs to maintain default factories that apply to multiple CVs.
  CompositeVectorSpace* Update(const CompositeVectorSpace& other);

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
                                     Teuchos::RCP<const Epetra_BlockMap> mastermap,
                                     Teuchos::RCP<const Epetra_BlockMap> ghostmap,
                                     int num_dofs);

  CompositeVectorSpace*
  AddComponents(const std::vector<std::string>& names,
                const std::vector<AmanziMesh::Entity_kind>& locations,
                std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> mastermaps,
                std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> ghostmaps,
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
                                     Teuchos::RCP<const Epetra_BlockMap> mastermap,
                                     Teuchos::RCP<const Epetra_BlockMap> ghostmap,
                                     int num_dof);

  CompositeVectorSpace*
  SetComponents(const std::vector<std::string>& names,
                const std::vector<AmanziMesh::Entity_kind> locations,
                std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> mastermaps,
                std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> ghostmaps,
                const std::vector<int>& num_dofs);

  // Write components to outstream.
  void Print(std::ostream& os) const
  {
    for (int i = 0; i != names_.size(); ++i) {
      os << "comp=" << names_[i] << " location=" << locations_[i] << " dofs=" << num_dofs_[i]
         << std::endl;
    }
  }

 private:
  // Indexing of name->int
  int Index_(const std::string& name) const
  {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
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
                           std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>>& mastermaps1,
                           std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>>& ghostmaps1,
                           std::vector<std::string>& names2,
                           std::vector<AmanziMesh::Entity_kind>& locations2,
                           std::vector<int>& num_dofs2,
                           std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>>& mastermaps2,
                           std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>>& ghostmaps2);

 private:
  bool ghosted_;

  // data
  bool owned_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::vector<std::string> names_;
  std::map<std::string, int> indexmap_;

  std::vector<AmanziMesh::Entity_kind> locations_;

  std::vector<int> num_dofs_;
  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> mastermaps_;
  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> ghostmaps_;

  friend class CompositeVector;
};

} // namespace Amanzi

#endif
