/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* **************************************************************** */
void TestLinearPressure(bool saturated) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  std::string wrm = (saturated) ? "saturated" : "van Genuchten";
  std::cout << "\nTest: Tensor Richards, a cube model: wrm=" << wrm << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  std::string xmlFileName = "test/flow_richards_tensor.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  if (saturated) 
    plist->sublist("PKs").sublist("flow")
        .sublist("water retention models").sublist("WRM for All")
        .set<std::string>("water retention model", "saturated");

  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  RCP<const AmanziMesh::Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Richards_PK* RPK = new Richards_PK(plist, "flow", S, soln);

  RPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // create Richards problem
  RPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();

  /* calculate the constant Darcy mass velocity */
  double rho = S->Get<double>("const_fluid_density");
  double mu = (*S->Get<CompositeVector>("viscosity_liquid").ViewComponent("cell"))[0][0];
  auto g = S->Get<AmanziGeometry::Point>("gravity");

  std::string passwd("flow");
  auto& perm = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");

  Point K(perm[0][0], perm[1][0], perm[2][0]);  // model the permeability tensor
  Point u0(1.0, 1.0, 1.0);
  Point v0(3);

  for (int i = 0; i < 3; i++) v0[i] = -u0[i] / K[i];
  v0 *= mu / rho;
  v0 += g * rho;
  std::cout << "rho=" << rho << "  mu=" << mu << std::endl;
  std::cout << "K=" << K << "  gravity=" << g << std::endl;
  std::cout << "grad(p)=" << v0 << std::endl;

  // solver the problem
  TI_Specs ti_specs;
  ti_specs.T0 = 0.0;
  ti_specs.dT0 = 1.0;
  ti_specs.T1 = 100.0;
  ti_specs.max_itrs = 400;

  AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  RPK->CommitStep(0.0, 1.0, S);  // dummy times for flow

  /* check accuracy */
  const auto& pressure = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
  const auto& flux = *S->Get<CompositeVector>("darcy_flux").ViewComponent("face");

  double err_p = 0.0, err_u = 0.0;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double p_exact = v0 * xc;
    // std::cout << c << " p_num=" << pressure[0][c] << " p_ex=" << p_exact << std::endl;
    err_p += pow(pressure[0][c] - p_exact, 2.0);
  }
  err_p = sqrt(err_p);

  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    const Point normal = mesh->face_normal(f);
    double f_exact = u0 * normal / rho;
    err_u += pow(flux[0][f] - f_exact, 2.0);
    // std::cout << f << " " << xf << "  flux_num=" << flux[0][f] << " f_ex=" << f_exact << std::endl;
  }
  err_u = sqrt(err_u);

  printf("Errors err_p=%8.3g  err_u=%8.3g\n", err_p, err_u);
  CHECK(err_p < 1e-8);
  CHECK(err_u < 1e-8);

  delete RPK;
}

TEST(FLOW_RICHARDS_ACCURACY) {
  TestLinearPressure(false);
  TestLinearPressure(true);
}
