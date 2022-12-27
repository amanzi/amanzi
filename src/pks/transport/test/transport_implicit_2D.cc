/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

  Dynamics is due to molecular diffusion, so the second-order
  and first-order schemes are identical.
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
#include "IO.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Amanzi::Transport
#include "TransportImplicit_PK.hh"

Epetra_MultiVector
runTest(int order, const std::string& linsolver)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTEST: implicit advance, spatial order=" << order << std::endl;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read and modify parameter list
  std::string xmlFileName("test/transport_implicit_2D.xml");
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  plist->sublist("PKs").sublist("transport").set<int>("spatial discretization order", order);
  plist->sublist("PKs")
    .sublist("transport")
    .sublist("time integrator")
    .set<std::string>("linear solver", linsolver);

  // create a mesh
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 6.0, 2.0, 192, 2);

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree =
    plist->sublist("cycle driver").sublist("pk_tree").sublist("transport");

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  TransportImplicit_PK TPK(pk_tree, plist, S, soln);

  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  // initialize a transport process kernel
  TPK.Initialize();

  // advance the state
  double t_old(0.0), t_new(0.0), dt(2.0e+3);
  while (t_new < 1.0e+5) {
    t_new = t_old + dt;

    if (comm->MyPID() == 0) std::cout << "\nCycle: T=" << t_old << " dT=" << dt << std::endl;
    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
  }

  // compute error in the final solution
  auto tcc = *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");

  double err(0.0);
  for (int c = 0; c < tcc.MyLength(); ++c) {
    const auto& xc = mesh->cell_centroid(c);
    err += fabs(tcc[0][c] - std::erfc(xc[0] / 2));
    // std::cout << xc << " " << tcc[0][c] << std::endl;
  }

  double err_tmp(err);
  mesh->get_comm()->SumAll(&err_tmp, &err, 1);
  err /= tcc.GlobalLength();
  std::cout << "Mean error in fracture: " << err << std::endl;
  CHECK(err < 1e-3);

  WriteStateStatistics(*S);

  return tcc;
}


TEST(IMPLICIT_TRANSPORT_2D)
{
  auto tcc1 = runTest(1, "PCG");
  auto tcc2 = runTest(2, "");

  tcc1.Update(1.0, tcc2, -1.0);
  double norm;
  tcc1.Norm2(&norm);
  CHECK_CLOSE(norm, 0.0, 2e-12);
}
