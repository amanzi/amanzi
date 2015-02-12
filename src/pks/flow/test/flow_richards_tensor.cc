/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "State.hh"
#include "Richards_PK.hh"


/* **************************************************************** */
TEST(FLOW_RICHARDS_ACCURACY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

std::cout << "Test: Tensor Richards, a cube model" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/flow_richards_tensor.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  RCP<const AmanziMesh::Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, gm);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  ParameterList state_list = plist.get<Teuchos::ParameterList>("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<Teuchos::ParameterList> global_list(&plist, Teuchos::RCP_WEAK_NO_DEALLOC);
  Richards_PK* RPK = new Richards_PK(global_list, "Flow", S);
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  RPK->InitializeFields();
  S->CheckAllFieldsInitialized();

  /* create Richards problem */
  RPK->Initialize(S.ptr());
  RPK->ti_specs_sss().T1 = 100.0;
  RPK->ti_specs_sss().max_itrs = 400;

  RPK->InitSteadyState(0.0, 1.0);

  /* calculate the constant Darcy mass velocity */
  double rho = *S->GetScalarData("fluid_density");
  double mu = *S->GetScalarData("fluid_viscosity");
  const AmanziGeometry::Point& g = RPK->gravity();

  std::string passwd("state");
  Epetra_MultiVector& perm = *S->GetFieldData("permeability", passwd)->ViewComponent("cell");

  Point K(perm[0][0], perm[1][0], perm[2][0]);  // model the permeability tensor
  Point u0(1.0, 1.0, 1.0);
  Point v0(3);

  for (int i = 0; i < 3; i++) v0[i] = -u0[i] / K[i];
  v0 *= mu / rho;
  v0 += g * rho;
  std::cout << "rho=" << rho << "  mu=" << mu << std::endl;
  std::cout << "K=" << K << "  gravity=" << g << std::endl;
  std::cout << "grad(p)=" << v0 << std::endl;

  RPK->AdvanceToSteadyState(0.0, 1.0);
  RPK->CommitState(0.0, S.ptr());

  /* check accuracy */
  const Epetra_MultiVector& pressure = *S->GetFieldData("pressure", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux", passwd)->ViewComponent("face");

  double err_p = 0.0, err_u = 0.0;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double p_exact = v0 * xc;
    // std::cout << c << " p_num=" << pressure[0][c] << " p_ex=" << p_exact << std::endl;
    err_p += pow(pressure[0][c] - p_exact, 2.0);
  }
  err_p = sqrt(err_p);

  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) {
    const Point& xf = mesh->face_centroid(f);
    const Point normal = mesh->face_normal(f);

    double p_exact = v0 * xf;
    double f_exact = u0 * normal / rho;
    err_u += pow(flux[0][f] - f_exact, 2.0);
    // std::cout << f << " " << xf << "  flux_num=" << flux[0][f] << " f_ex=" << f_exact << std::endl;
  }
  err_u = sqrt(err_u);

  printf("Errors err_p=%8.3g  err_u=%8.3g\n", err_p, err_u);
  CHECK(err_p < 1e-8);
  CHECK(err_u < 1e-8);

  delete comm;
  delete RPK;
}
