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
#include "MeshAudit.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"

/* **************************************************************** */
void runTest(double switch_time, std::string xmlfile, std::string exofile,
             std::string limiter, std::string stencil) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

std::cout << "Test: Advance on a 2D square mesh: limiter=" << limiter 
          << ", stencil=" << stencil << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  // read parameter list
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlfile);

  /* create a mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  RCP<const Mesh> mesh;
  if (exofile != "")
    mesh = meshfactory.create(exofile);
  else
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);

  plist->sublist("PKs").sublist("transport").sublist("reconstruction")
        .set<std::string>("limiter", limiter)
        .set<std::string>("limiter stencil", stencil);
  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup(S.ptr());
  TPK.CreateDefaultState(mesh, 2);
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  std::string passwd("state"); 
  Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(1.0, 1.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux[0][f] = velocity * normal;
  }

  // initialize a transport process kernel from a transport state
  TPK.Initialize(S.ptr());

  // advance the transport state 
  int iter;
  double t_old(0.0), t_new(0.0), dt;
  bool flag(true);

  Teuchos::RCP<Epetra_MultiVector> 
      tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  iter = 0;
  while (t_new < 0.5) {
    if (t_new > switch_time && flag) {
      flag = false;
      flux.Scale(-1.0);
      std::cout << "Changing Darcy velocity direction to opposite.\n\n";
    }

    dt = TPK.StableTimeStep();
    t_new = t_old + dt;

    S->set_initial_time(t_old);
    S->set_intermediate_time(t_old);
    S->set_final_time(t_new);

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, S);

    t_old = t_new;
    iter++;

    if (iter < 15) {
      printf("T=%8.4f  C_0(x):", t_new);
      for (int k = 0; k < 9; k++) {
        int k1 = 9 - k;  // reflects cell numbering in the exodus file
        printf("%7.4f", (*tcc)[0][k1]); 
      }
      printf("\n");
    }

    if (t_new < 0.15) {
      GMV::open_data_file(*mesh, (std::string)"transport.gmv");
      GMV::start_data();
      GMV::write_cell_data(*tcc, 0, "Component_0");
      GMV::write_cell_data(*tcc, 1, "Component_1");
      GMV::close_data_file();
    }
  }


  // check that the final state is constant for no swicth time
  if (switch_time > 0.5) {
    for (int k = 0; k < 10; k++) {
      CHECK_CLOSE(1.0, (*tcc)[0][k], 1e-6);
    }
  } else {
    for (int k = 0; k < 10; k++) {
      CHECK_CLOSE(0.0, (*tcc)[0][k], 2e-6);
    }
  }

  
}


TEST(ADVANCE_2D_MESH) {
  // no velocity switch
  runTest(1.0, "test/transport_2D.xml", "", "tensorial", "face to cells");
  runTest(1.0, "test/transport_2D.xml", "test/median7x8.exo", "Barth-Jespersen", "cell to closest cells");
}

TEST(ADVANCE_2D_MESH_SWITCH_FLOW) {
  runTest(0.16, "test/transport_2D.xml", "", "tensorial", "face to cells");
  runTest(0.16, "test/transport_2D.xml", "", "tensorial", "cell to closest cells");
}

TEST(ADVANCE_2D_MESH_SWITCH_FLOW_KUZMIN) {
  runTest(0.16, "test/transport_2D.xml", "", "Kuzmin", "node to cells");
}

