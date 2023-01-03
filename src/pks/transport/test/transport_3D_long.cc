/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport

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
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"

/* **************************************************************** */
TEST(ADVANCE_WITH_3D_MESH)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: 2.5D transport on a cubic mesh for long time" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_3D_long.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create a mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));
  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  pref.push_back(Framework::STK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create("test/rect3D_50x50x1.exo");

  //Amanzi::MeshAudit audit(mesh);
  //audit.Verify();

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = false;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);

  /* modify the default state for the problem at hand */
  std::string passwd("state");
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    flux[0][f] = velocity * normal;
  }

  /* initialize a transport process kernel */
  TPK.Initialize();

  /* advance the transport state */
  int iter;
  double t_old(0.0), t_new(0.0), dt;

  auto& tcc =
    *S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  iter = 0;
  bool flag = true;
  while (t_new < 0.3) {
    dt = TPK.StableTimeStep(-1);
    t_new = t_old + dt;

    S->set_initial_time(t_old);
    S->set_intermediate_time(t_old);
    S->set_final_time(t_new);

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    if (t_new > 0.1 && flag) {
      flag = false;
      if (TPK.MyPID == 0) {
        GMV::open_data_file(*mesh, (std::string) "transport.gmv");
        GMV::start_data();
        GMV::write_cell_data(tcc, 0, "Component_0");
        GMV::close_data_file();
      }
      break;
    }
  }
}
