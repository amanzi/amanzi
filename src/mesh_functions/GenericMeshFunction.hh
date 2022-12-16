/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Mesh Functions

  It is like a unique mesh function but returns one string value per domain.
*/

#ifndef AMANZI_GENERIC_MESH_FUNCTION_HH_
#define AMANZI_GENERIC_MESH_FUNCTION_HH_


#include <utility>
#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Mesh.hh"
#include "TabularStringFunction.hh"

namespace Amanzi {
namespace Functions {

template <class Generic>
class GenericMeshFunction {
 public:
  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;

  typedef std::pair<Teuchos::RCP<Domain>, Teuchos::RCP<const Generic>> Spec;
  typedef std::vector<Teuchos::RCP<Spec>> SpecList;

  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;
  typedef std::pair<Techous::RCP<Spec>, Teuchos::RCP<MeshIDs>> MeshSpec;
  typedef std::vector<Teuchos::RCP<MeshSpec>> MeshSpecList;

 public:
  GenericMeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh){};
  ~GenericMeshFunction(){};

  // method that check uniqueness
  void AddSpec(const Teuchos::RCP<Spec>& spec);

 protected:
  std::map<AmanziMesh::Entity_kind, Teuchos::RCP<MeshSpecList>> mesh_specs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

} // namespace Functions
} // namespace Amanzi


#endif
