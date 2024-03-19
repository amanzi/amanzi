/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include "stdlib.h"
#include "math.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "bilinear_form_reg.hh"
#include "CycleDriver.hh"
#include "IO.hh"
#include "MeshAudit.hh"
#include "eos_reg.hh"
#include "evaluators_flow_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mechanics_reg.hh"
#include "pks_mpc_reg.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

double AdaptiveIntegral(double (*f)(double, double*),
                        double *params,
                        double a, double b, double tol = 1e-6, int n = 10)
{
  double val1(0.0), val2(1.0e+98), h((b - a) / n);

  while (fabs(val2 - val1) > tol * fabs(val2 + 1e-10)) {
    val1 = val2;
    val2 = 0.0;
    for (int i = 0; i < n; ++i) {
      double x = a + i * h + h / 2; 
      val2 += h * f(x, params);
    }
    n *= 2;
    h /= 2;
  }

  return val2;
}


/*
  The McNmee-Gibson problme with analytic solution
*/

double fun(double u, double *params) {
  double t = params[0];
  double x = params[1];
  double y = params[2];
  double K = 2 / M_PI / u * std::cos(x * u) * std::sin(u);
  double t2 = std::sqrt(t);
  return K * (std::exp(-u * y) * (1 + std::erf(u * t2)) - exp(-u * y) * std::erfc(y / 2 / t2 - u * t2));
}


void
RunTest(const std::string xmlInFileName)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  int nx(30), ny(45);
  double L(6.0);
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, L, 9.0, nx, ny);

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::ObservationData obs_data;
  Amanzi::CycleDriver cd(plist, S, comm, obs_data);

  double dt(0.02);

  cd.Init_PK(0);
  cd.Setup();

  S->set_time(0.0);
  S->set_cycle(0);
  S->set_position(0);

  // compute initial pressure
  double biot(1.0), P0(7000.0), a(1.0);
  double rho(1e+3), mu(1e-3), k(1e-12), phi_cp(0.0);
  double factor(P0 / M_PI); 
  double params[3];

  auto& p_c = *S->GetW<CompositeVector>("pressure", Tags::DEFAULT, "").ViewComponent("cell");
  int ncells = p_c.MyLength();
  for (int c = 0; c < ncells; ++c) {
    const AmanziGeometry::Point xc = mesh->getCellCentroid(c);

    double x1(xc[0] - a), x2(xc[0] + a), y(9.0 - xc[1]);
    double tmp = std::asin(2 * a * y / std::sqrt((x1 * x1 + y * y) * (x2 * x2 + y * y)));
    if (x1 * x2 + y * y > 0.0) p_c[0][c] = factor * tmp;
    else p_c[0][c] = factor * (M_PI - tmp);
  }

  auto& p_f = *S->GetW<CompositeVector>("pressure", Tags::DEFAULT, "").ViewComponent("face");
  int nfaces = p_f.MyLength();
  for (int f = 0; f < nfaces; ++f) {
    const AmanziGeometry::Point xf = mesh->getFaceCentroid(f);

    double x1(xf[0] - a), x2(xf[0] + a), y(9.0 - xf[1]);
    double tmp = std::asin(2 * a * y / std::sqrt((x1 * x1 + y * y) * (x2 * x2 + y * y)));
    if (x1 * x2 + y * y > 0.0) p_f[0][f] = factor * tmp;
    else p_f[0][f] = factor * (M_PI - tmp);
  }
  Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(S->GetEvaluatorPtr("pressure", Tags::DEFAULT))->SetChanged();

  // initialization of PK
  cd.Initialize();

  cd.set_dt(dt);
  dt = cd.get_dt();

  S->GetW<double>("dt", Tags::NEXT, "dt") = dt;
  S->GetRecordW("dt", Tags::NEXT, "dt").set_initialized();

  S->InitializeEvaluators();
  S->CheckAllFieldsInitialized();
  S->InitializeIOFlags();

  cd.Visualize();
  cd.Observations();
  WriteStateStatistics(*S);

  // iterate time stepping
  while (S->get_time() < 2.0) {
    Utils::Units units("molar");
    std::cout << "\nCYCLE " << S->get_cycle()
              << ": time = " << units.OutputTime(S->get_time())
              << ", dt = " << units.OutputTime(dt) << "\n";

    S->GetW<double>("dt", Tags::DEFAULT, "dt") = dt;
    S->set_initial_time(S->get_time());
    S->set_final_time(S->get_time() + dt);
    S->set_intermediate_time(S->get_time());
    S->set_position(0);

    dt = cd.Advance(dt);
  }

  // compute analytic solution
  double perr(0.0), pnorm(0.0);
  std::string label = obs_data.observationLabels()[0];

  for (auto& quad : obs_data[label]) {
    params[0] = quad.time;
    params[1] = 0.1;
    params[2] = 2.1;
    double p = P0 * AdaptiveIntegral(&fun, params, 1.0e-10, 10.0, 1e-5, 100);
    // std::cout << quad.time << " p=" << p << " " << quad.value << std::endl;

    perr += std::pow(p - quad.value, 2.0);
    pnorm += std::pow(p, 2.0);
  }

  double tmp(perr);
  mesh->getComm()->SumAll(&tmp, &perr, 1);
  tmp = pnorm;
  mesh->getComm()->SumAll(&tmp, &pnorm, 1);

  perr = std::sqrt(perr / pnorm);
std::cout << perr << std::endl;
  CHECK(perr < 0.03);
}


TEST(MPC_DRIVER_BIOT_MCNAMEE_GIBSON)
{
  RunTest("test/mpc_biot_mcnamee_gibson.xml");
}

