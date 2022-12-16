/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! A factory for embedding Sperry-style root meshes within a standard 3D mesh.
/*!

This factory looks to embed Sperry-style root meshes within a 3D standard
domain subsurface mesh.  It optionally builds a logical mesh (for sequential
coupling) or an embedded mesh (for global implicit coupling).

We note that this class requires two steps:

1. Create() a topologically relevant mesh, and
2. SetGeometry() to specify the geometry of this mesh.

This is important as it opens the potential for doing dynamic plant traits.

*/

#ifndef AMANZI_LOGICAL_MESH_SPERRY_HH_
#define AMANZI_LOGICAL_MESH_SPERRY_HH_

#include "Mesh.hh"
#include "MeshLogicalFactory.hh"


namespace Amanzi {
namespace AmanziMesh {

class MeshLogicalEmbeddedSperry : public MeshLogicalFactory {
 public:
  MeshLogicalEmbeddedSperry(const Epetra_MpiComm* comm,
                      const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm) :
      MeshLogicalFactory(comm, gm)    
  {}

  Teuchos::RCP<MeshLogical> CreateLogical(int n_leaf, int n_stem,
          double max_rooting_depth, int n_rheizosphere_shells);
  Teuchos::RCP<MeshLogical> CreateEmbedded(int n_leaf, int n_stem,
          double max_rooting_depth, int n_rheizosphere_shells,
          const Teuchos::RCP<AmanziMesh::Mesh>& bg_mesh);

 protected:
  std::string Name_(const std::string& pftname, int col,
                    const std::string& component, int col_cell=0) {
    std::string name;
    if (component == "aroot" || component == "rheizosphere") {
      return pftname+"_"+std::to_string(col)+"_"+std::to_string(col_cell)+"_"+component;
    } else {
      return pftname+"_"+std::to_string(col)+"_"+component;
    }
  }
  
 protected:
  // topological parameters
  int n_leaf_, n_stem_;
  double max_rooting_depth_;
  int n_shells_;
  Teuchos::RCP<AmanziMesh::Mesh> bg_mesh_;

  std::vector<Entity_ID_List> rheizosphere_to_bg_;
  
};



} // namespace AmanziMesh
} // namespace Amanzi

#endif
