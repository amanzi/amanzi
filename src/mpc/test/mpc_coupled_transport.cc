/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <utility>
#include "math.h"
#include "stdlib.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "eos_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"

std::pair<Epetra_MultiVector, Epetra_MultiVector>
RunTest(int order, double cfl)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nTEST: coupled transport, implicit scheme order=" << order << std::endl;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::string xmlInFileName = "test/mpc_coupled_transport.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  plist->sublist("PKs").sublist("transport matrix").set<int>("spatial discretization order", order);
  plist->sublist("PKs")
    .sublist("transport fracture")
    .set<int>("spatial discretization order", order);

  plist->sublist("PKs")
    .sublist("transport matrix")
    .sublist("reconstruction")
    .set<double>("limiter cfl", cfl);
  plist->sublist("PKs")
    .sublist("transport fracture")
    .sublist("reconstruction")
    .set<double>("limiter cfl", cfl);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  int n(12);
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 0.0, 6.0, 6.0, 6.0, n, n, n);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture = factory.create(mesh, names, AmanziMesh::Entity_kind::FACE);

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  // verify bounds and 1D-symmetry
  const auto& tcc_m =
    *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  double cmin(1e+99), cmax(-1e+99), xc0(6.0 - 3.0 / n);
  for (int c = 0; c < n * n * n; ++c) {
    const auto& xc = mesh->getCellCentroid(c);
    if (std::fabs(xc[0] - xc0) < 1e-3) {
      cmin = std::min(cmin, tcc_m[0][c]);
      cmax = std::max(cmin, tcc_m[0][c]);
    }
  }
  CHECK_CLOSE(cmin, cmax, 1e-12);

  const auto& tcc_f =
    *S->Get<CompositeVector>("fracture-total_component_concentration").ViewComponent("cell");
  cmin = 1e+99;
  cmax = -1e+99;
  for (int c = 0; c < n * n; ++c) {
    CHECK(tcc_f[0][c] <= 1.0);
    CHECK(tcc_f[0][c] >= 0.0);
    const auto& xc = mesh_fracture->getCellCentroid(c);
    if (std::fabs(xc[0] - xc0) < 1e-3) {
      cmin = std::min(cmin, tcc_f[0][c]);
      cmax = std::max(cmin, tcc_f[0][c]);
    }
  }
  CHECK_CLOSE(cmin, cmax, 1e-12);

  return std::make_pair(tcc_f, tcc_m);
}


TEST(MPC_DRIVER_COUPLED_TRANSPORT)
{
  auto tcc1 = RunTest(1, 1.0);
  auto tcc2 = RunTest(2, 1.0);

  double err, norm;
  tcc1.first.Update(1.0, tcc2.first, -1.0);
  tcc1.first.Norm2(&err);
  std::cout << "Error between 1st and 2nd: " << err << std::endl;

  tcc1.second.Norm2(&norm);
  tcc1.second.Update(1.0, tcc2.second, -1.0);
  tcc1.second.Norm2(&err);
  std::cout << "Error between 1st and 2nd: " << err << " norm: " << norm << std::endl;
  CHECK_CLOSE(err / norm, 0.15, 0.1);
}
