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


namespace ATS {
namespace Testing {

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical>
plantMesh(const Teuchos::RCP<Epetra_MpiComm>& comm,
          const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
          bool include_soil);

// Teuchos::RCP<Amanzi::AmanziMesh::MeshEmbeddedLogical>
// plantMeshEmbedded();


} // namespace
} // namespace

#endif
