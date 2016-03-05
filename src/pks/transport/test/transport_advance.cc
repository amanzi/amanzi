/*
  The transport component of the Amanzi code, serial unit tests.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "FrameworkTraits.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "Transport_PK.hh"


TEST(ADVANCE_WITH_MESH_FRAMEWORK) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::vector<std::string> framework_name;
  framework_name.push_back("MSTK");
  framework_name.push_back("STK");
  framework_name.push_back("SIMPLE");
  Framework framework[3] = {MSTK, STKMESH, Simple};   

  for (int frm = 0; frm < 3; frm++) {
    std::cout << "Test: advance with framework " << framework_name[frm] << std::endl;
    if(!framework_available(framework[frm])) continue;
#ifdef HAVE_MPI
    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

    // read parameter list
    std::string xmlFileName("test/transport_advance.xml");
    if (frm == 2) xmlFileName = "test/transport_advance_simple.xml";
    Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

    // create a mesh
    ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, comm));

    FrameworkPreference pref;
    pref.clear();
    pref.push_back(framework[frm]);

    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
    RCP<const Mesh> mesh;
    if (frm < 2) {
      mesh = meshfactory("test/hex_3x3x3_ss.exo", gm);
    } else {
      mesh = meshfactory(0.0,0.0,0.0, 1.0,1.0,1.0, 20, 1, 1, gm); 
    }
  
    // create a simple state and populate it
    Amanzi::VerboseObject::hide_line_prefix = false;

    std::vector<std::string> component_names;
    component_names.push_back("Component 0");
    component_names.push_back("Component 1");

    Teuchos::ParameterList state_list = plist->sublist("State");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    S->set_intermediate_time(0.0);

    Transport_PK TPK(plist, S, "Transport", component_names);
    TPK.Setup(S.ptr());
    TPK.CreateDefaultState(mesh, 2);
    S->InitializeFields();
    S->InitializeEvaluators();

    // modify the default state for the problem at hand
    std::string passwd("state"); 
    Teuchos::RCP<Epetra_MultiVector> 
        flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

    AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f = 0; f < nfaces_owned; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      (*flux)[0][f] = velocity * normal;
    }

    // initialize a transport process kernel
    TPK.Initialize(S.ptr());

    // advance the state
    double t_old(0.0), t_new, dt;
    dt = TPK.CalculateTransportDt();
    t_new = t_old + dt;
    TPK.AdvanceStep(t_old, t_new);

    // printing cell concentration
    Teuchos::RCP<Epetra_MultiVector> 
        tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell");

    while(t_new < 1.2) {
      dt = TPK.CalculateTransportDt();
      t_new = t_old + dt;

      TPK.AdvanceStep(t_old, t_new);
      TPK.CommitStep(t_old, t_new, S);

      t_old = t_new;
 
      if (t_new < 0.4) {
        printf("T=%6.2f  C_0(x):", t_new);
        for (int k = 0; k < 9; k++) printf("%7.4f", (*tcc)[0][k]); std::cout << std::endl;
      }
    }

    // check that the final state is constant
    for (int k = 0; k < 4; k++) 
      CHECK_CLOSE((*tcc)[0][k], 1.0, 1e-6);

    if (frm == 2) {
      for (int k = 0; k < 19; k++) {
         CHECK(((*tcc)[0][k] - (*tcc)[0][k+1]) > -1e-15);
      }
    }

    delete comm;
  }
}
 


