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
#include "Point.hh"

#include "State.hh"
#include "Transport_PK.hh"


double f_step(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if (x[0] <= 1 + t) return 1;
  return 0;
}

double f_smooth(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

double f_cubic(const Amanzi::AmanziGeometry::Point& x, double t) {
  if( x[0] < 1 + t ) return 1;
  if( x[0] > 3 + t ) return 0;
  double z = (x[0]-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}


/* **************************************************************** */
TEST(DISPERSION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: dispersion" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/transport_dispersion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create an MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);
  int nx = 20;
  RCP<const Mesh> mesh = factory(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 10, 1, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Transport_PK TPK(plist, S, component_names);
  TPK.CreateDefaultState(mesh, 1);

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Teuchos::RCP<Epetra_MultiVector> 
      flux = S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    (*flux)[0][f] = velocity * normal;
  }

  Teuchos::RCP<Epetra_MultiVector> 
      tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*tcc)[0][c] = f_step(xc, 0.0);
  }

  S->GetFieldData("porosity", passwd)->PutScalar(1.0);
  *(S->GetScalarData("fluid_density", passwd)) = 1.0;

  /* initialize a transport process kernel */
  Amanzi::VerboseObject::hide_line_prefix = true;
  TPK.InitPK();
  TPK.PrintStatistics();

  /* advance the state */
  double dT, dT0;
  dT0 = TPK.CalculateTransportDt();

  int i, k, iter = 0;
  double T = 0.0, T1 = 1.0;
  while (T < T1) {
    dT = std::min(TPK.CalculateTransportDt(), T1 - T);
    dT = std::min(dT, dT0);
    // for (int k = 0; k < nx; k++) printf("%10.8f\n", (*tcc)[0][k]); 
    // printf("\n");

    TPK.Advance(dT);
    TPK.CommitState(S);
    T += dT;
    iter++;

    TPK.CheckTracerBounds(*tcc, 0, 0.0, 1.0, 1e-12);
  }
 
  delete comm;
}



