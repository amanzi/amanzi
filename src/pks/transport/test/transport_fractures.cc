/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "State.hh"

// Transport
#include "Transport_PK.hh"

/* **************************************************************** */
TEST(ADVANCE_TWO_FRACTURES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

std::cout << "Test: Advance on a 2D square mesh" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  // read parameter list
  std::string xmlFileName = "test/transport_fractures.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, comm));

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({Framework::MSTK, Framework::STKMESH}));
  RCP<const Mesh> mesh3D = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, gm);

  // extract fractures mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  RCP<const Mesh_MSTK> mesh_mstk = rcp_static_cast<const Mesh_MSTK>(mesh3D);
  RCP<const Mesh>  mesh = Teuchos::rcp(new Mesh_MSTK(*mesh_mstk, setnames, AmanziMesh::FACE));

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // create a simple state and populate it
  Amanzi::VerboseObject::hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  Transport_PK TPK(plist, S, "transport", component_names);
  TPK.Setup(S.ptr());
  TPK.CreateDefaultState(mesh, 1);
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state
  Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux_fracture", "state")->ViewComponent("cell");

  int dir;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  AmanziGeometry::Point velocity(1.0, 0.2, -0.1);
  for (int c = 0; c < ncells_owned; c++) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh->face_normal(f, false, c, &dir);
      flux[n][c] = velocity * normal;
    }
  }
  S->GetField("darcy_flux_fracture", "state")->set_initialized();

  // we still need the old flux until testing is complete
  Epetra_MultiVector& flux_old = *S->GetFieldData("darcy_flux", "state")->ViewComponent("face");
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux_old[0][f] = velocity * normal;
  }

  // initialize the transport process kernel
  TPK.Initialize(S.ptr());

  // advance the transport state 
  int iter, k;
  double t_old(0.0), t_new(0.0), dt;
  Epetra_MultiVector& tcc = *S->GetFieldData("total_component_concentration", "state")
                              ->ViewComponent("cell", false);

  iter = 0;
  while (t_new < 0.2) {
    dt = TPK.StableTimeStep();
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, S);

    t_old = t_new;
    iter++;

    // verify solution
    for (int c = 0; c < ncells_owned; ++c) 
        CHECK(tcc[0][c] >= 0.0 && tcc[0][c] <= 1.0);
  }

  // test the maximum principle
  AmanziMesh::Entity_ID_List block;
  mesh->get_set_entities("fracture 2", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

  // test that solute enter the second fracture
  double tcc_max(0.0);
  for (int n = 0; n < block.size(); ++n) {
    tcc_max = std::max(tcc_max, tcc[0][block[n]]);
  }
  CHECK(tcc_max > 0.25);

  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(tcc, 0, "Component_0");
  GMV::close_data_file();

  delete comm;
}





