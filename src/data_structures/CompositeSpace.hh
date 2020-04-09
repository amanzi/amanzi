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
  A CompositeSpace is a BlockSpace where all maps come from a specific mesh.

  This should be thought of as a vector-space: it lays out data components as a
  mesh along with entities on the mesh.  This meta-data can be used with the
  mesh's *_map() methods to create the data.

  This class is very light weight as it maintains only meta-data.
*/

#ifndef AMANZI_COMPOSITE_SPACE_HH_
#define AMANZI_COMPOSITE_SPACE_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "BlockSpace.hh"

namespace Amanzi {

namespace AmanziMesh {
class Mesh;
}


// Nonmember helper function
std::pair<Map_ptr_type, Map_ptr_type>
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location);


class CompositeSpace : public BlockSpace {
 public:
  // Constructor
  CompositeSpace(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                 const std::vector<std::string>& names,
                 const std::map<std::string, AmanziMesh::Entity_kind> locations,
                 const std::map<std::string, BlockMap_ptr_type>& master_maps,
                 const std::map<std::string, BlockMap_ptr_type>& ghost_maps,
                 const std::map<std::string, std::size_t>& num_vectors,
                 bool ghosted);

  CompositeSpace(const CompositeSpace& other) = default;
  CompositeSpace& operator=(const CompositeSpace&) = default;

  //
  // Specs for the construction of CVs
  // -------------------------------------
  // mesh specification
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return mesh_; }
  AmanziMesh::Entity_kind Location(const std::string& name) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::map<std::string, AmanziMesh::Entity_kind> locations_;
  bool ghosted_;

  friend class CompositeVectorSpace;
};

} // namespace Amanzi

#endif
