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
#include "Mesh_simple.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "Richards_PK.hpp"
#include "BDF2_Dae.hpp"


/* **************************************************************** */
TEST(FLOW_RICHARDS_ACCURACY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

cout << "Test: Tensor Richards, a cube model" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_tensor.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  RCP<AmanziMesh::Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,0.0, 1.0,1.0,1.0, 2, 2, 2, comm, gm)); 

  // create the state
  ParameterList state_list = parameter_list.get<Teuchos::ParameterList>("State");
  State* S = new State(state_list, mesh);
  RCP<Flow_State> FS = Teuchos::rcp(new AmanziFlow::Flow_State(*S));

  // create Richards problem
  ParameterList flow_list = parameter_list.get<Teuchos::ParameterList>("Flow");
  Richards_PK* RPK = new Richards_PK(flow_list, FS);
  RPK->set_standalone_mode(true);
  RPK->InitPK();
  RPK->InitSteadyState(0.0, 1e-8);

  // calculate the constant Darcy mass velocity
  double rho = FS->ref_fluid_density();
  double mu = FS->ref_fluid_viscosity();
  AmanziGeometry::Point& g = RPK->gravity();

  const Epetra_Vector& kh = FS->ref_horizontal_permeability();
  const Epetra_Vector& kv = FS->ref_vertical_permeability();

  Point K(kh[0], kh[0], kv[0]);  // model the permeability tensor
  Point u0(1.0, 1.0, 1.0);
  Point v0(3);

  for (int i=0; i<3; i++) v0[i] = -u0[i] / K[i];
  v0 *= mu / rho;
  v0 += g * rho;
  cout << "rho=" << rho << "  mu=" << mu << endl;
  cout << "K=" << K << "  gravity=" << g << endl;
  cout << "grad(p)=" << v0 << endl;

  RPK->AdvanceToSteadyState();
  RPK->CommitStateForTransport(FS);

  // check accuracy
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
 
  double err_p = 0.0, err_u = 0.0;
  int ncells = mesh->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double p_exact = v0 * xc;
    cout << c << " p_num=" << pressure[c] << " p_ex=" << p_exact << endl;
    err_p += pow(pressure[c] - p_exact, 2.0);
  }
  err_p = sqrt(err_p);

  int nfaces = mesh->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f<nfaces; f++) {
    const Point& xf = mesh->face_centroid(f);
    const Point normal = mesh->face_normal(f);
  
    double p_exact = v0 * xf;
    double f_exact = u0 * normal / rho;
    err_u += pow(darcy_flux[f] - f_exact, 2.0);
    //cout << f << " " << xf << "  flux_num=" << darcy_flux[f] << " f_ex=" << f_exact << endl;
  }
  err_u = sqrt(err_u);

  CHECK(err_p < 1e-8);
  CHECK(err_u < 1e-8);
 
  delete comm;
  delete RPK;
  delete S;
}
