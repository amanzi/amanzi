/*
The flow component of the Amanzi code, richards unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_MSTK.hh"
#include "MeshAudit.hh"
#include "gmv_mesh.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "Darcy_PK.hpp"


/* **************************************************************** */
TEST(FLOW_2D_DARCY_WELL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "Test: 2D specific storage Darcy, homogeneous media" << endl;

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_darcy_well.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);
  RCP<Mesh> mesh = rcp(new Mesh_MSTK(-10.0,-5.0, 10.0,0.0, 200,50, &comm, gm));

  // create and populate flow state
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(mesh));
  FS->set_permeability(0.1, 2.0, "Computational domain");
  FS->set_porosity(0.2);
  FS->set_fluid_viscosity(1.0);
  FS->set_fluid_density(1.0);
  FS->set_gravity(-1.0);

  // create Richards process kernel
  ParameterList flow_list = parameter_list.get<ParameterList>("Flow");
  Darcy_PK* DPK = new Darcy_PK(flow_list, FS);
  DPK->set_standalone_mode(true);
  DPK->InitPK();
  DPK->InitSteadyState(0.0, 1e-8);

  // create the initial pressure function
  Epetra_Vector& p = FS->ref_pressure();

  for (int c=0; c<p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[c] = xc[1] * (xc[1] + 2.0);
  }

  // transient solution
  double dT = 1.0;
  for (int n=0; n<50; n++) {
    DPK->advance(dT); 
    DPK->commitState(FS);
 
    if (MyPID == 0) {
      GMV::open_data_file(*mesh, (std::string)"flow.gmv");
      GMV::start_data();
      GMV::write_cell_data(FS->ref_pressure(), 0, "pressure");
      GMV::close_data_file();
    }
  }

  delete DPK;
}
