/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
      Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

/*
  Factory for a CompositeSpace, which is then used to create CompositeVectors.
*/

#include <algorithm>

#include "dbc.hh"
#include "errors.hh"

#include "Mesh.hh"
#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"

namespace Amanzi {


// constructor
CompositeVectorSpace::CompositeVectorSpace()
  : owned_(false), mesh_(Teuchos::null), ghosted_(false){};

// constructor from Space
CompositeVectorSpace::CompositeVectorSpace(const CompositeSpace& map)
  : owned_(false), mesh_(map.Mesh()), ghosted_(map.ghosted_)
{
  for (const auto& name : map) {
    names_.emplace_back(name);
    locations_.emplace_back(map.Location(name));
    num_vectors_.emplace_back(map.getNumVectors(name));
    mastermaps_[name] = map.ComponentMap(name, false);
    ghostmaps_[name] = map.ComponentMap(name, true);
  }

  // idiot check!
  AMANZI_ASSERT(CreateSpace()->SameAs(map));
}


// Check equivalence of spaces.
bool
CompositeVectorSpace::SameAs(const CompositeVectorSpace& other) const
{
  auto cvs = this->CreateSpace();
  auto other_cvs = other.CreateSpace();
  if (cvs != Teuchos::null && other_cvs != Teuchos::null)
    return cvs->SameAs(*other_cvs);
  return SubsetOf(other) && other.SubsetOf(*this);
}

// Check equivalence of spaces.
bool
CompositeVectorSpace::LocallySameAs(const CompositeVectorSpace& other) const
{
  auto cvs = this->CreateSpace();
  auto other_cvs = other.CreateSpace();
  if (cvs != Teuchos::null && other_cvs != Teuchos::null)
    return cvs->LocallySameAs(*other_cvs);
  return SubsetOf(other) && other.SubsetOf(*this);
}

// Check subset of spaces.
bool
CompositeVectorSpace::SubsetOf(const CompositeVectorSpace& other) const
{
  if (mesh_ != other.mesh_) return false;
  for (const auto& name : *this) {
    if (!other.HasComponent(name)) return false;
    if (getNumVectors(name) != other.getNumVectors(name)) return false;
    if (Location(name) != other.Location(name)) return false;
  }
  return true;
}


// Create just the CompositeSpace
Teuchos::RCP<const CompositeSpace>
CompositeVectorSpace::CreateSpace() const
{
  if (mesh_ == Teuchos::null) {
    Errors::Message message(
      "CompositeVectorSpace: Cannot create a Space prior to setting a Mesh.");
    throw(message);
  }
  if (size() == 0) {
    Errors::Message message(
      "CompositeVectorSpace: Cowardly refusing to make an empty map.");
    throw(message);
  }

  // TODO: really need to refactor and update the internal data, but for now
  // just reorganize. --etc
  std::map<std::string, AmanziMesh::Entity_kind> map_locs;
  std::map<std::string, std::size_t> map_num_vecs;
  for (std::size_t i = 0; i != size(); ++i) {
    map_locs[names_[i]] = locations_[i];
    map_num_vecs[names_[i]] = num_vectors_[i];
  }

  // Create the map
  return Teuchos::rcp(new CompositeSpace(
    mesh_, names_, map_locs, mastermaps_, ghostmaps_, map_num_vecs, ghosted_));
}


Comm_ptr_type
CompositeVectorSpace::Comm() const
{
  if (mesh_ == Teuchos::null)
    return Teuchos::null;
  else
    return mesh_->get_comm();
}

// -------------------------------------
// Specs for the construction of CVs
// -------------------------------------

// Update all specs from another factory's specs.
// Useful for PKs to maintain default factories that apply to multiple CVs.
CompositeVectorSpace*
CompositeVectorSpace::Update(const CompositeVectorSpace& other)
{
  if (this != &other) {
    AddComponents(other.names_,
                  other.locations_,
                  other.mastermaps_,
                  other.ghostmaps_,
                  other.num_vectors_);
    SetMesh(other.mesh_);
  }
  return this;
};


// ghosted
CompositeVectorSpace*
CompositeVectorSpace::SetGhosted(bool ghosted)
{
  ghosted_ |= ghosted;
  return this;
}

// component specification
CompositeVectorSpace*
CompositeVectorSpace::SetMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  if (mesh_ == Teuchos::null) {
    mesh_ = mesh;
  } else if (mesh_ != mesh) {
    Errors::Message message("CompositeVectorSpace: SetMesh called on space "
                            "that already has a mesh with a different mesh.");
    throw(message);
  }
  return this;
};


// Add methods append their specs to the factory's spec, checking to make
// sure the spec is OK if the full spec has been set (by an owning PK).
CompositeVectorSpace*
CompositeVectorSpace::AddComponent(const std::string& name,
                                   AmanziMesh::Entity_kind location,
                                   int num_dof)
{
  std::vector<std::string> names(1, name);
  std::vector<AmanziMesh::Entity_kind> locations(1, location);
  std::vector<int> num_dofs(1, num_dof);
  AddComponents(names, locations, num_dofs);
  return this;
};

CompositeVectorSpace*
CompositeVectorSpace::AddComponents(
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::vector<int>& num_dofs)
{
  std::map<std::string, BlockMap_ptr_type> mastermaps;
  std::map<std::string, BlockMap_ptr_type> ghostmaps;

  for (int i = 0; i < locations.size(); ++i) {
    auto maps = getMaps(*mesh_, locations[i]);
    mastermaps[names[i]] = maps.first;
    ghostmaps[names[i]] = maps.second;
  }

  return AddComponents(names, locations, mastermaps, ghostmaps, num_dofs);
}

CompositeVectorSpace*
CompositeVectorSpace::AddComponent(const std::string& name,
                                   AmanziMesh::Entity_kind location,
                                   const BlockMap_ptr_type& mastermap,
                                   const BlockMap_ptr_type& ghostmap,
                                   int num_dof)
{
  std::vector<std::string> names(1, name);
  std::vector<int> num_dofs(1, num_dof);
  std::vector<AmanziMesh::Entity_kind> locations(1, location);
  std::map<std::string, BlockMap_ptr_type> mastermaps;
  std::map<std::string, BlockMap_ptr_type> ghostmaps;

  mastermaps[name] = mastermap;
  ghostmaps[name] = ghostmap;

  return AddComponents(names, locations, mastermaps, ghostmaps, num_dofs);
}


CompositeVectorSpace*
CompositeVectorSpace::AddComponents(
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::map<std::string, BlockMap_ptr_type>& mastermaps,
  const std::map<std::string, BlockMap_ptr_type>& ghostmaps,
  const std::vector<int>& num_dofs)
{
  // Add a set of specifications to the factory's list.
  if (owned_) {
    // Factory's specs are fixed by an owning PK.  Check that the requested
    // components are supplied by the owning PK.
    if (!CheckContained_(names_, names)) {
      Errors::Message message;
      message << "Requested components {";
      for (const auto& n : names) {
        message << "\"" << n << "\",";
      }
      message << "} are not supplied by an already owned CompositeVector"
              << " which has components {";
      for (const auto& n : names_) {
        message << "\"" << n << "\",";
      }
      message << "}.";
      Exceptions::amanzi_throw(message);
    }

  } else {
    // Factory's specs are not yet fixed.  Form the union (checking
    // consistency) and save this as the factory's new specs.
    if (!UnionAndConsistent_(names,
                             locations,
                             num_dofs,
                             mastermaps,
                             ghostmaps,
                             names_,
                             locations_,
                             num_vectors_,
                             mastermaps_,
                             ghostmaps_)) {
      Errors::Message message(
        "Requested components are not consistent with previous request.");
      Exceptions::amanzi_throw(message);
    }
  }

  InitIndexMap_();

  return this;
};


// Set methods fix the component specs, checking to make sure all previously
// added specs are contained in the new spec.
CompositeVectorSpace*
CompositeVectorSpace::SetComponent(const std::string& name,
                                   AmanziMesh::Entity_kind location,
                                   int num_dof)
{
  std::vector<std::string> names(1, name);
  std::vector<AmanziMesh::Entity_kind> locations(1, location);
  std::vector<int> num_dofs(1, num_dof);
  SetComponents(names, locations, num_dofs);
  return this;
};


CompositeVectorSpace*
CompositeVectorSpace::SetComponents(
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::vector<int>& num_dofs)
{
  std::map<std::string, BlockMap_ptr_type> mastermaps;
  std::map<std::string, BlockMap_ptr_type> ghostmaps;

  for (int i = 0; i < locations.size(); ++i) {
    mastermaps[names[i]] = mesh_->map(locations[i], false);
    ghostmaps[names[i]] = mesh_->map(locations[i], true);
  }

  return SetComponents(names, locations, mastermaps, ghostmaps, num_dofs);
};

CompositeVectorSpace*
CompositeVectorSpace::SetComponent(const std::string& name,
                                   AmanziMesh::Entity_kind location,
                                   const BlockMap_ptr_type& mastermap,
                                   const BlockMap_ptr_type& ghostmap,
                                   int num_dof)
{
  std::vector<std::string> names(1, name);
  std::vector<int> num_dofs(1, num_dof);
  std::vector<AmanziMesh::Entity_kind> locations(1, location);
  std::map<std::string, BlockMap_ptr_type> mastermaps;
  std::map<std::string, BlockMap_ptr_type> ghostmaps;

  mastermaps[name] = mastermap;
  ghostmaps[name] = ghostmap;

  return SetComponents(names, locations, mastermaps, ghostmaps, num_dofs);
}

CompositeVectorSpace*
CompositeVectorSpace::SetComponents(
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::map<std::string, BlockMap_ptr_type>& mastermaps,
  const std::map<std::string, BlockMap_ptr_type>& ghostmaps,
  const std::vector<int>& num_dofs)
{
  if (owned_) {
    // check equal
    if (names != names_ || num_dofs != num_vectors_) {
      Errors::Message message(
        "SetComponent() called on an owned space with a differing spec.");
      Exceptions::amanzi_throw(message);
    }
    return this;
  }

  // Make sure everything we've been asked for can be covered by this spec.
  if (!CheckContained_(names, names_)) {
    Errors::Message message(
      "CompositeVector's components do not supply a previous request.");
    Exceptions::amanzi_throw(message);
  }

  // Make sure the specs are consistent.
  if (!CheckConsistent_(
        names_, locations_, num_vectors_, names, locations, num_dofs)) {
    Errors::Message message("CompositeVector's components are not consistent "
                            "with a previous request.");
    Exceptions::amanzi_throw(message);
  }

  // Re-spec to the new spec and declare ourselves owned.
  names_ = names;
  locations_ = locations;
  num_vectors_ = num_dofs;
  mastermaps_ = mastermaps;
  ghostmaps_ = ghostmaps;

  owned_ = true;
  InitIndexMap_();

  return this;
};


// Set up the map from names to index into vector.
void
CompositeVectorSpace::InitIndexMap_()
{
  for (int i = 0; i != names_.size(); ++i) { indexmap_[names_[i]] = i; }
}

BlockMap_ptr_type
CompositeVectorSpace::getMap(const std::string& name, bool ghost) const
{
  if (std::find(names_.begin(), names_.end(), name) == names_.end()) {
    Errors::Message message("Map: Requested component (" + name +
                            ") doesn't exist in CompositeVectorSpace.");
    Exceptions::amanzi_throw(message);
  }

  if (ghost) {
    return ghostmaps_.at(name);
  } else {
    return mastermaps_.at(name);
  }
}


// Assert that one spec is contained in another.
bool
CompositeVectorSpace::CheckContained_(
  const std::vector<std::string>& containing,
  const std::vector<std::string>& contained)
{
  for (auto name = contained.begin(); name != contained.end(); ++name) {
    if (std::find(containing.begin(), containing.end(), *name) ==
        containing.end()) {
      // special case if name is "boundary_face", then check for "face" and use
      // Vandelay
      if (*name == std::string("boundary_face")) {
        if (std::find(containing.begin(),
                      containing.end(),
                      std::string("face")) == containing.end()) {
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
bool
CompositeVectorSpace::CheckConsistent_(
  const std::vector<std::string>& names1,
  const std::vector<AmanziMesh::Entity_kind>& locations1,
  const std::vector<int>& num_dofs1, const std::vector<std::string>& names2,
  const std::vector<AmanziMesh::Entity_kind>& locations2,
  const std::vector<int>& num_dofs2)
{
  for (int i = 0; i != names1.size(); ++i) {
    std::vector<std::string>::const_iterator n2_it =
      std::find(names2.begin(), names2.end(), names1[i]);
    // special case for Vandelay
    if (n2_it == names2.end() && names1[i] == std::string("boundary_face")) {
      n2_it = std::find(names2.begin(), names2.end(), std::string("face"));
    }

    if (n2_it != names2.end()) {
      int j = n2_it - names2.begin();
      if (locations1[i] != locations2[j]) { return false; }
      if (num_dofs1[i] != num_dofs2[j]) { return false; }
    }
  }
  return true;
};


// Form the union, returning it in the 2 vectors.
bool
CompositeVectorSpace::UnionAndConsistent_(
  const std::vector<std::string>& names1,
  const std::vector<AmanziMesh::Entity_kind>& locations1,
  const std::vector<int>& num_dofs1, std::vector<std::string>& names2,
  std::vector<AmanziMesh::Entity_kind>& locations2, std::vector<int>& num_dofs2)
{
  for (int i = 0; i != names1.size(); ++i) {
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
          if (num_dofs1[i] != num_dofs2[j]) { return false; }
        } else {
          // add this spec
          names2.push_back(names1[i]);
          locations2.push_back(locations1[i]);
          num_dofs2.push_back(num_dofs1[i]);
        }
      } else if (names1[i] == std::string("face")) {
        // check if names1 is face, but names2 has boundary_face
        n2_it =
          std::find(names2.begin(), names2.end(), std::string("boundary_face"));
        if (n2_it != names2.end()) {
          int j = n2_it - names2.begin();
          if ((locations1[i] != AmanziMesh::FACE) ||
              (locations2[j] != AmanziMesh::BOUNDARY_FACE)) {
            return false;
          }
          if (num_dofs1[i] != num_dofs2[j]) { return false; }

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
      if (locations1[i] != locations2[j]) { return false; }
      if (num_dofs1[i] != num_dofs2[j]) { return false; }
    }
  }
  return true;
};


// Form the union, returning it in the 2 vectors.
bool
CompositeVectorSpace::UnionAndConsistent_(
  const std::vector<std::string>& names1,
  const std::vector<AmanziMesh::Entity_kind>& locations1,
  const std::vector<int>& num_dofs1,
  const std::map<std::string, BlockMap_ptr_type>& mastermaps1,
  const std::map<std::string, BlockMap_ptr_type>& ghostmaps1,
  std::vector<std::string>& names2,
  std::vector<AmanziMesh::Entity_kind>& locations2, std::vector<int>& num_dofs2,
  std::map<std::string, BlockMap_ptr_type>& mastermaps2,
  std::map<std::string, BlockMap_ptr_type>& ghostmaps2)
{
  for (int i = 0; i != names1.size(); ++i) {
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
          if (num_dofs1[i] != num_dofs2[j]) { return false; }
        } else {
          // add this spec
          names2.push_back(names1[i]);
          locations2.push_back(locations1[i]);
          num_dofs2.push_back(num_dofs1[i]);
          mastermaps2[names1[i]] = mastermaps1.at(names1[i]);
          ghostmaps2[names1[i]] = ghostmaps1.at(names1[i]);
        }
      } else if (names1[i] == std::string("face")) {
        // check if names1 is face, but names2 has boundary_face
        n2_it =
          std::find(names2.begin(), names2.end(), std::string("boundary_face"));
        if (n2_it != names2.end()) {
          int j = n2_it - names2.begin();
          if ((locations1[i] != AmanziMesh::FACE) ||
              (locations2[j] != AmanziMesh::BOUNDARY_FACE)) {
            return false;
          }
          if (num_dofs1[i] != num_dofs2[j]) { return false; }

          // union of face and boundary_face is face
          names2[j] = "face";
          locations2[j] = AmanziMesh::FACE;
        } else {
          // add this spec
          names2.push_back(names1[i]);
          locations2.push_back(locations1[i]);
          num_dofs2.push_back(num_dofs1[i]);
          mastermaps2[names1[i]] = mastermaps1.at(names1[i]);
          ghostmaps2[names1[i]] = ghostmaps1.at(names1[i]);
        }
      } else {
        // add this spec
        names2.push_back(names1[i]);
        locations2.push_back(locations1[i]);
        num_dofs2.push_back(num_dofs1[i]);
        mastermaps2[names1[i]] = mastermaps1.at(names1[i]);
        ghostmaps2[names1[i]] = ghostmaps1.at(names1[i]);
      }

    } else {
      // just make sure they match
      int j = n2_it - names2.begin();
      if (locations1[i] != locations2[j]) { return false; }
      if (num_dofs1[i] != num_dofs2[j]) { return false; }
    }
  }
  return true;
};


} // namespace Amanzi
