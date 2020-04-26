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
#include "UniqueLocalIndex.hh"

// Transport
#include "TransportExplicit_PK.hh"

/* **************************************************************** */
TEST(ADVANCE_TWO_FRACTURES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

std::cout << "Test: Advance on a 2D square mesh" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  // read parameter list
  std::string xmlFileName = "test/transport_fractures.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  RCP<const Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  // extract fractures mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  // RCP<const Mesh> mesh = meshfactory.create(mesh3D, setnames, AmanziMesh::FACE);
  RCP<const Mesh> mesh = meshfactory.create("test/fractures.exo");

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  std::cout << "pid=" << comm->MyPID() << " cells: " << ncells_owned 
                                       << " faces: " << nfaces_owned << std::endl;

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(0.0);

  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup(S.ptr());
  TPK.CreateDefaultState(mesh, 1);
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state
  auto& flux = *S->GetFieldData("darcy_flux", "state")->ViewComponent("face", true);
  const auto flux_map = S->GetFieldData("darcy_flux", "state")->Map().Map("face", true);

  int dir;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  AmanziGeometry::Point velocity(1.0, 0.2, -0.1);
  for (int c = 0; c < ncells_wghost; c++) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      int g = flux_map->FirstPointInElement(f);
      int ndofs = flux_map->ElementSize(f);
      if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh, f, c);

      const AmanziGeometry::Point& normal = mesh->face_normal(f, false, c, &dir);
      flux[0][g] = (velocity * normal) * dir;
    }
  }
  S->GetField("darcy_flux", "state")->set_initialized();

  // initialize the transport process kernel
  TPK.Initialize(S.ptr());

  // advance the transport state 
  int iter;
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
  double tmp = tcc_max;
  mesh->get_comm()->MaxAll(&tmp, &tcc_max, 1);
  CHECK(tcc_max > 0.25);

  GMV::open_data_file(*mesh, (std::string)"transport.gmv");
  GMV::start_data();
  GMV::write_cell_data(tcc, 0, "Component_0");
  GMV::close_data_file();
}





