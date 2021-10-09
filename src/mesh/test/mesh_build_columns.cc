/*
  Copyright 2010- held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

/*
  build_columns() insantiates extra data structures in a Mesh that provide
  semi-structured functionality, specifically around meshes whose vertical
  structure consists of stacks of cells, or columns.  Fully unstructured in x
  and y, nodes and cells are vertically structured.

  Note that pinchouts and variances in dz cause cell centroids to not be
  vertically stacked!  But nodes will always be vertically stacked.  In fact,
  this can be extreme enough that a cell's centroid may have a larger z
  coordinate than the face above its centroid.

  This test works with an extruded mesh, where the dz is uniform across the
  layer of cells, in which case cell centroids are also vertically stacked.

*/

#include <vector>
#include <iostream>

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshAudit.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

using namespace Amanzi;

TEST(MESH_COLUMNS)
{
  auto comm = getDefaultComm();
  int nprocs = comm->NumProc();

  // We are not including MOAB since Mesh_MOAB.cc does not have
  // routines for generating a mesh
  std::vector<AmanziMesh::Framework> frameworks;
  if (AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
  }
  if (nprocs == 1) {
    frameworks.push_back(AmanziMesh::Framework::SIMPLE);
  }

  auto mesh_pars = Teuchos::rcp(new Teuchos::ParameterList());
  mesh_pars->set<std::string>("partitioner", "zoltan_rcb");

  for (const auto& frm : frameworks) {
    std::cerr << "Testing columns with " << AmanziMesh::to_string(frm) << std::endl;

    // Create the mesh
    auto mesh = createStructuredUnitHex({frm}, 3,3,10, comm, Teuchos::null, mesh_pars, 2,3,2);

    // Explicitly call build columns method
    mesh->buildColumns();

    // call test function
    testColumnsUniformDz(*mesh, 2.0/10);
  }
}

