//
// Set of functions that create demo meshes for testing.
//
// ------------------------------------------------------------------

#ifndef MESH_LOGICAL_DEMO_MESHES_HH_
#define MESH_LOGICAL_DEMO_MESHES_HH_

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshEmbeddedLogical.hh"

namespace Amanzi {
namespace Testing {

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalSegmentRegularManual();

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalSegmentIrregularManual();

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalYManual();

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalY();

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
demoMeshLogicalFromXML(const std::string& meshname);

Teuchos::RCP<Amanzi::AmanziMesh::MeshEmbeddedLogical>
demoMeshLogicalYEmbedded();
  
class RegularMeshCellFromCoordFunctor {
 public:
  RegularMeshCellFromCoordFunctor(const Amanzi::AmanziGeometry::Point& X0,
				  const Amanzi::AmanziGeometry::Point& X1,
				  int nx, int ny, int nz);

  Amanzi::AmanziMesh::Entity_ID operator()(const Amanzi::AmanziGeometry::Point& p);

 private:
  Amanzi::AmanziGeometry::Point X0_;
  Amanzi::AmanziGeometry::Point X1_;
  Amanzi::AmanziGeometry::Point dX_;
  int nx_,ny_,nz_;
};

} // namespace Testing
} // namespace Amanzi

#endif
