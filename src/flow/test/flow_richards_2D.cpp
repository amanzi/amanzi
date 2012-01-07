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

#include "MSTK_types.h"
#include "Mesh_MSTK.hh"
#include "MeshAudit.hh"
#include "gmv_mesh.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "Richards_PK.hpp"
#include "BDF2_Dae.hpp"


/* **************************************************************** */
TEST(FLOW_2D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

cout << "Test: 2D Richards, 2-layer model" << endl;
  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_2D.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list);
  RCP<Mesh> mesh = rcp(new Mesh_MSTK(0.0,-2.0, 1.0,0.0, 10,40, MPI_COMM_WORLD, gm)); 

  // create flow state
  ParameterList& state_list = parameter_list.get<ParameterList>("State");
  State S(state_list, mesh);
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));

  // create Richards process kernel
  ParameterList flow_list = parameter_list.get<ParameterList>("Flow");
  ParameterList rp_list = flow_list.get<ParameterList>("Richards Problem");
  Richards_PK* RPK = new Richards_PK(rp_list, FS);
  RPK->Init();

  // create the initial pressure function
  Epetra_Vector p(RPK->get_super_map());
  Epetra_Vector* pcells = FS->create_cell_view(p);
  Epetra_Vector* pfaces = FS->create_face_view(p);

  for (int c=0; c<pcells->MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    (*pcells)[c] = xc[1] * (xc[1] + 2.0);
  }

  for (int f=0; f<pfaces->MyLength(); f++) {
    const Point& xf = mesh->face_centroid(f);
    (*pfaces)[f] = xf[1] * (xf[1] + 2.0);
  }

  S.update_pressure(*pcells);
  S.set_time(0.0);

  // solve the problem
  RPK->advance_to_steady_state();
 
  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(RPK->get_solution_cells(), 0, "pressure");
  GMV::close_data_file();

  delete RPK;
}
