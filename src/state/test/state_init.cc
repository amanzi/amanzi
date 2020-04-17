/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "CompositeVector.hh"
#include "MeshFactory.hh"

// State
#include "State.hh"

TEST(FIELD_INITIALIZATION)
{
  using namespace Amanzi;

  // Epetra_MpiComm comm(MPI_COMM_WORLD);
  auto comm = getDefaultComm();
  std::string xmlFileName = "test/state_init.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create geometric model
  Teuchos::ParameterList region_list =
    plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create a mesh
  AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto mesh = meshfactory.create("test/cube3x3x3.exo");

  // create a state
  Teuchos::ParameterList state_list =
    plist.get<Teuchos::ParameterList>("state");
  State S(state_list);

  // populate state
  S.RegisterMesh("domain", mesh);

  // constant value
  S.Require<CompositeVector, CompositeVectorSpace>("field1")
    .SetMesh(mesh)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  //
  S.Require<CompositeVector, CompositeVectorSpace>("field2")
    .SetMesh(mesh)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>("permeability")
    .SetMesh(mesh)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 3);

  // setup creates data
  S.Setup();

  // initialize data
  S.InitializeFields();

  // check state's fields
  // -- porosity (simple field)
  int ncells =
    mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto phi =
    S.Get<CompositeVector>("field1").ViewComponent<MirrorHost>("cell");
  for (int c = 0; c < ncells; ++c) { CHECK_EQUAL(0.25, phi(c, 0)); }

  // from exo currently not supported in new state
  // // -- scalar field from a file
  // const Epetra_MultiVector& fld =
  // *S.Get<CompositeVector>("field").ViewComponent("cell"); for (int c = 1; c <
  // ncells; ++c) {
  //   CHECK_EQUAL(1.0e-11, fld[0][c]);
  // }
  // // -- permeability from a file
  // const Epetra_MultiVector& K =
  // *S.Get<CompositeVector>("permeability").ViewComponent("cell"); for (int c =
  // 1; c < ncells; ++c) {
  //   CHECK_EQUAL(1.0e-11, K[0][c]);
  //   CHECK_EQUAL(1.0e-12, K[1][c]);
  //   CHECK_EQUAL(1.0e-13, K[2][c]);
  // }
}
