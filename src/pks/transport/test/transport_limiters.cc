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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"

#include "State.hh"
#include "Transport_PK.hh"
#include "TransportExplicit_PK.hh"


/* ****************************************************************
 * Test LimiterBarthJespersen()routine.
 * ************************************************************* */
TEST(LIMITER_BARTH_JESPERSEN)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: read transport XML file" << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_limiters.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory factory(comm);
  factory.set_preference(pref);
  RCP<const Mesh> mesh = factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 4, 7, gm);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  S->set_intermediate_time(0.0);


  TransportExplicit_PK TPK(plist, S, component_names);
  TPK.CreateDefaultState(mesh, 1);

  /* modify the default state for the problem at hand */
  std::string passwd("state");
  auto flux = S->GetFieldData("volumetric_flow_rate", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  /* initialize a transport process kernel */
  TPK.Initialize(S.ptr());
  double dT = TPK.CalculateTransportDt(); // We call it to identify upwind cells.

  /* create a linear field */
  const Epetra_Map& cmap = mesh->cell_map(false);
  RCP<Epetra_Vector> scalar_field = rcp(new Epetra_Vector(cmap));

  CompositeVectorSpace cv_space;
  cv_space.SetMesh(mesh);
  cv_space.SetGhosted(true);
  cv_space.SetComponent("cell", AmanziMesh::CELL, 3);

  RCP<CompositeVector> gradient =
    Teuchos::RCP<CompositeVector>(new CompositeVector(cv_space, true));
  RCP<Epetra_MultiVector> grad = gradient->ViewComponent("cell", false);

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*scalar_field)[c] = 5.0 - xc[0] - 0.5 * xc[1] - 0.2 * xc[2];
    (*grad)[0][c] = -1.0;
    (*grad)[1][c] = -0.5;
    (*grad)[2][c] = -0.2;
  }

  /* calculate and verify limiters */
  RCP<Epetra_Vector> limiter = rcp(new Epetra_Vector(cmap));
  TPK.LimiterBarthJespersen(0, scalar_field, gradient, limiter);

  for (int c = 0; c < ncells - 1; c++) { // the corner cell gives limiter=0
    CHECK_CLOSE(1.0, (*limiter)[c], 1e-6);
  }
}
