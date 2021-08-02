/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <vector>
#include <iostream>

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "Geometry.hh"
#include "Mesh.hh"
#include "mesh_mstk/Mesh_MSTK.hh"
#include "MeshFactory.hh"
#include "MeshException.hh"

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
  mesh_pars->sublist("unstructured").sublist("expert").set<std::string>("partitioner", "zoltan_rcb");

  for (const auto& frm : frameworks) {
    std::cerr << "Testing columns with " << AmanziMesh::framework_names.at(frm) << std::endl;

    // Create the mesh
    auto mesh = createFrameworkStructuredUnitHex({frm}, 3,3,10, comm, Teuchos::null, mesh_pars, 2,3,2);

    // Explicitly call build columns method
    mesh->build_columns();

    // call test function
    testColumnsUniformDz(mesh, 2.0/10);
  }
}

