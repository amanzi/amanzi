/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Factory for a CompositeVector on an Amanzi mesh.
   ------------------------------------------------------------------------- */
#include <algorithm>

#include "dbc.hh"
#include "errors.hh"

#include "CompositeVectorSpace.hh"

namespace Amanzi {

// constructor
CompositeVectorSpace::CompositeVectorSpace() :
  owned_(false), mesh_(Teuchos::null), ghosted_(false) {}

// copy constructor
CompositeVectorSpace::CompositeVectorSpace(const CompositeVectorSpace& other) :
    ghosted_(other.ghosted_),
    owned_(other.owned_),
    mesh_(other.mesh_),
    names_(other.names_),
    indexmap_(other.indexmap_),
    locations_(other.locations_),
    num_dofs_(other.num_dofs_) {}

// copy constructor
CompositeVectorSpace::CompositeVectorSpace(const CompositeVectorSpace& other,
        bool ghosted) :
    ghosted_(ghosted),
    owned_(other.owned_),
    mesh_(other.mesh_),
    names_(other.names_),
    indexmap_(other.indexmap_),
    locations_(other.locations_),
    num_dofs_(other.num_dofs_) {}


// Check equivalence of spaces.
bool
CompositeVectorSpace::SameAs(const CompositeVectorSpace& other) const {
  if (NumComponents() != other.NumComponents()) return false;
  if (mesh_ != other.mesh_) return false;
  for (name_iterator name=begin(); name!=end(); ++name) {
    if (!other.HasComponent(*name)) return false;
    if (NumVectors(*name) != other.NumVectors(*name)) return false;
    if (Location(*name) !=other.Location(*name)) return false;
  }
  return true;
}

// Check subset of spaces.
bool
CompositeVectorSpace::SubsetOf(const CompositeVectorSpace& other) const {
  if (mesh_ != other.mesh_) return false;
  for (name_iterator name=begin(); name!=end(); ++name) {
    if (!other.HasComponent(*name)) return false;
    if (NumVectors(*name) != other.NumVectors(*name)) return false;
    if (Location(*name) !=other.Location(*name)) return false;
  }
  return true;
}

// -------------------------------------
// Specs for the construction of CVs
// -------------------------------------

// Update all specs from another factory's specs.
// Useful for PKs to maintain default factories that apply to multiple CVs.
CompositeVectorSpace*
CompositeVectorSpace::Update(const CompositeVectorSpace& other) {
  if (other.owned_) {
    SetComponents(other.names_, other.locations_, other.num_dofs_);
  } else {
    AddComponents(other.names_, other.locations_, other.num_dofs_);
  }

  SetMesh(other.mesh_);
  //  SetGhosted(other.ghosted_);
  return this;
};


// ghosted
CompositeVectorSpace*
CompositeVectorSpace::SetGhosted(bool ghosted) {
  ghosted_ |= ghosted;
  return this;
}

// component specification
CompositeVectorSpace*
CompositeVectorSpace::SetMesh(
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  if (mesh_ == Teuchos::null) {
    mesh_ = mesh;
  } else {
    ASSERT(mesh_ == mesh);
  }
  return this;
};


// Add methods append their specs to the factory's spec, checking to make
// sure the spec is OK if the full spec has been set (by an owning PK).
CompositeVectorSpace*
CompositeVectorSpace::AddComponent(std::string name,
        AmanziMesh::Entity_kind location,
        int num_dof) {
  std::vector<std::string> names(1,name);
  std::vector<AmanziMesh::Entity_kind> locations(1,location);
  std::vector<int> num_dofs(1,num_dof);
  AddComponents(names, locations, num_dofs);
  return this;
};


CompositeVectorSpace*
CompositeVectorSpace::AddComponents(const std::vector<std::string>& names,
                   const std::vector<AmanziMesh::Entity_kind>& locations,
                   const std::vector<int>& num_dofs) {
  // Add a set of specifications to the factory's list.
  if (owned_) {
    // Factory's specs are fixed by an owning PK.  Check that the requested
    // components are supplied by the owning PK.
    if (!CheckContained_(names_, names)) {
      Errors::Message message("Requested components are not supplied by an already owned CompositeVector.");
      Exceptions::amanzi_throw(message);
    }

  } else {
    // Factory's specs are not yet fixed.  Form the union (checking
    // consistency) and save this as the factory's new specs.
    if (!UnionAndConsistent_(names, locations, num_dofs,
                             names_, locations_, num_dofs_)) {
      Errors::Message message("Requested components are not consistent with previous request.");
      Exceptions::amanzi_throw(message);
    }
  }

  InitIndexMap_();
  return this;
};


// Set methods fix the component specs, checking to make sure all previously
// added specs are contained in the new spec.
CompositeVectorSpace*
CompositeVectorSpace::SetComponent(std::string name,
        AmanziMesh::Entity_kind location,
        int num_dof) {
  std::vector<std::string> names(1,name);
  std::vector<AmanziMesh::Entity_kind> locations(1,location);
  std::vector<int> num_dofs(1,num_dof);
  SetComponents(names, locations, num_dofs);
  return this;
};


CompositeVectorSpace*
CompositeVectorSpace::SetComponents(const std::vector<std::string>& names,
        const std::vector<AmanziMesh::Entity_kind>& locations,
        const std::vector<int>& num_dofs) {
  // These components will be provided by an owning PK.
  if (owned_) {
    Errors::Message message("SetComponent() cannot be called by a non-owning PK, and this factory is already owned.");
    Exceptions::amanzi_throw(message);
  }

  // Make sure everything we've been asked for can be covered by this spec.
  if (!CheckContained_(names, names_)) {
    Errors::Message message("CompositeVector's components do not supply a previous request.");
    Exceptions::amanzi_throw(message);
  }

  // Make sure the specs are consistent.
  if (!CheckConsistent_(names_, locations_, num_dofs_,
                        names, locations, num_dofs)) {
      Errors::Message message("CompositeVector's components are not consistent with a previous request.");
      Exceptions::amanzi_throw(message);
  }

  // Re-spec to the new spec and declare ourselves owned.
  names_.clear();
  locations_.clear();
  num_dofs_.clear();

  names_ = names;
  locations_ = locations;
  num_dofs_ = num_dofs;

  owned_ = true;
  InitIndexMap_();
  return this;
};

// Set up the map from names to index into vector.
void CompositeVectorSpace::InitIndexMap_() {
  for (int i=0; i!=names_.size(); ++i) {
    indexmap_[names_[i]] = i;
  }
}


// Assert that one spec is contained in another.
bool CompositeVectorSpace::CheckContained_(const std::vector<std::string>& containing,
                     const std::vector<std::string>& contained) {
  for (std::vector<std::string>::const_iterator name=contained.begin();
       name!=contained.end(); ++name) {
    if (std::find(containing.begin(), containing.end(), *name) == containing.end()) {
      // special case if name is "boundary_face", then check for "face" and use Vandelay
      if (*name == std::string("boundary_face")) {
        if (std::find(containing.begin(), containing.end(), std::string("face")) == containing.end()) {
          return false;
        }
      } else {
        return false;
      }
    }
  }
  return true;
};


// Assert that the two specs are consistent.
bool CompositeVectorSpace::CheckConsistent_(const std::vector<std::string>& names1,
                      const std::vector<AmanziMesh::Entity_kind>& locations1,
                      const std::vector<int>& num_dofs1,
                      const std::vector<std::string>& names2,
                      const std::vector<AmanziMesh::Entity_kind>& locations2,
                      const std::vector<int>& num_dofs2) {
  for (int i=0; i!=names1.size(); ++i) {
    std::vector<std::string>::const_iterator n2_it =
        std::find(names2.begin(), names2.end(), names1[i]);
    // special case for Vandelay
    if (n2_it == names2.end() && names1[i] == std::string("boundary_face")) {
      n2_it = std::find(names2.begin(), names2.end(), std::string("face"));
    }

    if (n2_it != names2.end()) {
      int j = n2_it - names2.begin();
      if (locations1[i] != locations2[j]) {
        return false;
      }
      if (num_dofs1[i] != num_dofs2[j]) {
        return false;
      }
    }
  }
  return true;
};


// Form the union, returning it in the 2 vectors.
bool CompositeVectorSpace::UnionAndConsistent_(const std::vector<std::string>& names1,
                         const std::vector<AmanziMesh::Entity_kind>& locations1,
                         const std::vector<int>& num_dofs1,
                         std::vector<std::string>& names2,
                         std::vector<AmanziMesh::Entity_kind>& locations2,
                         std::vector<int>& num_dofs2) {
  for (int i=0; i!=names1.size(); ++i) {
    std::vector<std::string>::iterator n2_it =
      std::find(names2.begin(), names2.end(), names1[i]);
    if (n2_it == names2.end()) {
      // no match

      if (names1[i] == std::string("boundary_face")) {
        // check if names1 is boundary_face, but names2 has face
        n2_it = std::find(names2.begin(), names2.end(), std::string("face"));
        if (n2_it != names2.end()) {
          int j = n2_it - names2.begin();
          if ((locations1[i] != AmanziMesh::BOUNDARY_FACE) ||
              (locations2[j] != AmanziMesh::FACE)) {
            return false;
          }
          if (num_dofs1[i] != num_dofs2[j]) {
            return false;
          }
        } else {
          // add this spec
          names2.push_back(names1[i]);
          locations2.push_back(locations1[i]);
          num_dofs2.push_back(num_dofs1[i]);
        }
      } else if (names1[i] == std::string("face")) {
        // check if names1 is face, but names2 has boundary_face
        n2_it = std::find(names2.begin(), names2.end(), std::string("boundary_face"));
        if (n2_it != names2.end()) {
          int j = n2_it - names2.begin();
          if ((locations1[i] != AmanziMesh::FACE) ||
              (locations2[j] != AmanziMesh::BOUNDARY_FACE)) {
            return false;
          }
          if (num_dofs1[i] != num_dofs2[j]) {
            return false;
          }

          // union of face and boundary_face is face
          names2[j] = "face";
          locations2[j] = AmanziMesh::FACE;
        } else {
          // add this spec
          names2.push_back(names1[i]);
          locations2.push_back(locations1[i]);
          num_dofs2.push_back(num_dofs1[i]);
        }
      } else {
        // add this spec
        names2.push_back(names1[i]);
        locations2.push_back(locations1[i]);
        num_dofs2.push_back(num_dofs1[i]);
      }

    } else {
      // just make sure they match
      int j = n2_it - names2.begin();
      if (locations1[i] != locations2[j]) {
        return false;
      }
      if (num_dofs1[i] != num_dofs2[j]) {
        return false;
      }
    }
  }
  return true;
};


} // namespace Amanzi
