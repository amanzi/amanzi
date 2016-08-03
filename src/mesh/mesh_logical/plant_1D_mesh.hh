//
// Set of functions that create simple plant meshes.
//
// ------------------------------------------------------------------

#ifndef MESH_LOGICAL_PLANT_MESHES_HH_
#define MESH_LOGICAL_PLANT_MESHES_HH_

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshEmbeddedLogical.hh"


namespace Amanzi {
namespace Testing {

Teuchos::RCP<AmanziMesh::MeshLogical>
plantMesh(const Epetra_MpiComm* comm,
          const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
          bool include_soil);

// Teuchos::RCP<Amanzi::AmanziMesh::MeshEmbeddedLogical>
// plantMeshEmbedded();


} // namespace
} // namespace

#endif
