/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Factory for a CompositeVector on an Amanzi mesh.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITEVECTOR_FACTORY_HH_
#define AMANZI_COMPOSITEVECTOR_FACTORY_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "composite_vector.hh"

namespace Amanzi {

class CompositeVectorFactory {

public:
  // Constructor
  CompositeVectorFactory();
  CompositeVectorFactory(const CompositeVectorFactory& other);

  // Create the final product, a CompositeVector.
  Teuchos::RCP<CompositeVector> CreateVector(bool ghosted) const;

  // Is this factory and the resulting CV owned by a PK?
  bool owned() { return owned_; }
  void set_owned(bool owned=true) { owned_ = owned; }

  // -------------------------------------
  // Specs for the construction of CVs
  // -------------------------------------

  // Update all specs from another factory's specs.
  // Useful for PKs to maintain default factories that apply to multiple CVs.
  CompositeVectorFactory* Update(const CompositeVectorFactory& other);

  // ghost spec
  CompositeVectorFactory* SetGhosted(bool ghosted=true);
  bool Ghosted() { return ghosted_; }

  // mesh specification
  CompositeVectorFactory* SetMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return mesh_; }

  // component specification

  // Add methods append their specs to the factory's spec, checking to make
  // sure the spec is OK if the full spec has been set (by an owning PK).
  CompositeVectorFactory*
  AddComponent(std::string name,
               AmanziMesh::Entity_kind location,
               int num_dofs);

  CompositeVectorFactory*
  AddComponents(const std::vector<std::string>& names,
                const std::vector<AmanziMesh::Entity_kind>& locations,
                const std::vector<int>& num_dofs);

  // Set methods fix the component specs, checking to make sure all previously
  // added specs are contained in the new spec.
  CompositeVectorFactory*
  SetComponent(std::string name,
               AmanziMesh::Entity_kind location,
               int num_dofs);

  CompositeVectorFactory*
  SetComponents(const std::vector<std::string>& names,
                const std::vector<AmanziMesh::Entity_kind>& locations,
                const std::vector<int>& num_dofs);

private:
  // private and unimplemented
  CompositeVectorFactory& operator=(const CompositeVectorFactory&) {}

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

private:
  // data
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  bool owned_;

  // component information
  bool ghosted_;
  std::vector<std::string> names_;
  std::vector<AmanziMesh::Entity_kind> locations_;
  std::vector<int> num_dofs_;

};

} // namespace amanzi

#endif
