/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
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
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "eos_registration.hh"
#include "mdm_transport_registration.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "NumericalIntegration.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_transport_registration.hh"
#include "State.hh"


double
Case3(double a, double b, double c, double t)
{
  return a * pow(t, -1.5) * std::exp(-b / t + c * t);
}


double
RunTest(int icase,
        double u_f,
        double mol_diff_f,
        double mol_diff_m,
        double normal_diff,
        double L,
        double tend = 1e+5)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::string xmlInFileName = "test/mpc_coupled_dispersive_transport.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // modify the test
  plist->sublist("state")
    .sublist("initial conditions")
    .sublist("fracture-volumetric_flow_rate")
    .sublist("function")
    .sublist("All domain")
    .sublist("function")
    .sublist("dof 1 function")
    .sublist("function-constant")
    .set<double>("value", u_f);

  std::vector<double> tmp_f({ mol_diff_f });
  plist->sublist("PKs")
    .sublist("transport fracture")
    .sublist("molecular diffusion")
    .set<Teuchos::Array<double>>("aqueous values", tmp_f);

  std::vector<double> tmp_m({ mol_diff_m });
  plist->sublist("PKs")
    .sublist("transport matrix")
    .sublist("molecular diffusion")
    .set<Teuchos::Array<double>>("aqueous values", tmp_m);

  plist->sublist("state")
    .sublist("initial conditions")
    .sublist("fracture-normal_diffusion")
    .sublist("function")
    .sublist("All")
    .sublist("function")
    .sublist("function-constant")
    .set<double>("value", normal_diff);

  if (icase == 3) {
    plist->sublist("PKs")
      .sublist("transport fracture")
      .sublist("boundary conditions")
      .sublist("concentration")
      .sublist("tracer")
      .sublist("BC 0")
      .sublist("boundary concentration")
      .sublist("function-exprtk")
      .set<std::string>("formula", "exp(-1e-4 * t)");

    plist->sublist("cycle driver")
      .sublist("time periods")
      .sublist("TP 0")
      .set<double>("end period time", tend);
  }

  if (icase == 4) {
    plist->sublist("PKs")
      .sublist("transport fracture")
      .sublist("boundary conditions")
      .sublist("concentration")
      .sublist("tracer")
      .sublist("BC 0")
      .sublist("boundary concentration")
      .sublist("function-exprtk")
      .set<std::string>("formula", "if(t > 1.0, 0.0, 1.0)");

    plist->sublist("cycle driver")
      .sublist("time periods")
      .sublist("TP 0")
      .set<double>("end period time", tend);
  }

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 3.0 - L, 7.0, 1.0, 3.0 + L, 56, 1, 32, true, true);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture = factory.create(mesh, names, AmanziMesh::FACE);

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  const auto& tcc_m =
    *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  double cmin(1e+99), cmax(-1e+99);
  for (int c = 0; c < tcc_m.MyLength(); ++c) {
    const auto& xc = mesh->cell_centroid(c);
    if (std::fabs(xc[0] - 4.0625) < 1e-3) {
      cmin = std::min(cmin, tcc_m[0][c]);
      cmax = std::max(cmax, tcc_m[0][c]);
    }
    CHECK(tcc_m[0][c] < 1.0);
  }

  const auto& tcc_f =
    *S->Get<CompositeVector>("fracture-total_component_concentration").ViewComponent("cell");
  double fmin(1e+99), fmax(-1e+99), err(0.0), norm(0.0);

  double b = (*S->Get<CompositeVector>("fracture-aperture").ViewComponent("cell"))[0][0];
  double kn = (*S->Get<CompositeVector>("fracture-normal_diffusion").ViewComponent("cell"))[0][0];
  double phi = (*S->Get<CompositeVector>("porosity").ViewComponent("cell"))[0][0];

  double t(tend), Df(mol_diff_f), Dm(mol_diff_m);

  for (int c = 0; c < tcc_f.MyLength(); ++c) {
    const auto& xc = mesh_fracture->cell_centroid(c);
    if (std::fabs(xc[0] - 4.0625) < 1e-3) {
      fmin = std::min(fmin, tcc_f[0][c]);
      fmax = std::max(fmax, tcc_f[0][c]);
    }

    norm += fabs(tcc_f[0][c]);

    if (icase == 1) {
      double tmp2 = std::erfc(xc[0] / std::sqrt(Df * t) / 2);
      err += fabs(tcc_f[0][c] - tmp2);
      // std::cout << xc[0] << " " << tcc_f[0][c] << " " << tmp2 <<std::endl;
    } else if (icase == 2) {
      double Tm = u_f * t / b - xc[0];
      if (Tm > 0.0) {
        double tmp2 = std::erfc(xc[0] * kn * std::sqrt(b / (u_f * Dm * Tm)) / 2);
        err += fabs(tcc_f[0][c] - tmp2);
        // std::cout << xc[0] << " " << tcc_f[0][c] << " " << tmp2 <<std::endl;
      }
    } else if (icase == 3) {
      double Tm = t - xc[0] * b / u_f;
      if (Tm > 0.0) {
        double lambda(1e-4), a1, a2;
        a1 = xc[0] * phi * std::sqrt(Dm / M_PI) / u_f;
        a2 = a1 * a1 * M_PI;

        double c0 = std::exp(-lambda * Tm);
        double t1, h, tmp1, tmp2(0.0);
        for (int n = 1000; n < 1000000; n *= 10) {
          tmp1 = tmp2;
          tmp2 = 0.0;
          h = Tm / n;
          for (int i = 0; i < n; ++i) {
            t1 = i * h + h / 2;
            tmp2 += h * Case3(a1, a2, lambda, t1);
          }
          if (fabs(tmp1 - tmp2) < 1e-4 * (tmp1 + tmp2)) break;
        }
        err += fabs(tcc_f[0][c] - c0 * tmp2);
        // std::cout << xc[0] << " " << tcc_f[0][c] << " " << c0 * tmp2 <<std::endl;
      }
    } else {
      // do nothing
    }
  }
  double err_tmp(err), norm_tmp(norm);
  mesh->get_comm()->SumAll(&err_tmp, &err, 1);
  mesh->get_comm()->SumAll(&norm_tmp, &norm, 1);
  err /= tcc_f.GlobalLength();
  norm /= tcc_f.GlobalLength();
  std::cout << "Mean error in fracture: " << err << " solution norm=" << norm << std::endl;

  if (fmin < 1.0e+98) CHECK_CLOSE(fmin, fmax, 1e-12);


  // mass of solute in matrix
  // err = 0.0;
  norm = 0.0;
  for (int c = 0; c < tcc_m.MyLength(); ++c) {
    const auto& xc = mesh->cell_centroid(c);
    double vol = mesh->cell_volume(c);

    if (icase == 3) {
      double Tm = t - xc[0] * b / u_f;
      if (Tm > 0.0 && xc[2] > 3.0) {
        double lambda(1e-4), a1, a2;
        a1 = xc[0] * phi * std::sqrt(Dm / M_PI) / u_f + (xc[2] - 3.0 - b / 2) / sqrt(Dm * M_PI) / 2;
        a2 = a1 * a1 * M_PI;

        double c0 = std::exp(-lambda * Tm);
        double t1, h, tmp1, tmp2(0.0);
        for (int n = 1000; n < 1000000; n *= 10) {
          tmp1 = tmp2;
          tmp2 = 0.0;
          h = Tm / n;
          for (int i = 0; i < n; ++i) {
            t1 = i * h + h / 2;
            tmp2 += h * Case3(a1, a2, lambda, t1);
          }
          if (fabs(tmp1 - tmp2) < 1e-4 * (tmp1 + tmp2)) break;
        }
        // err += fabs(tcc_m[0][c] - c0 * tmp2);
        norm += c0 * tmp2 * vol;
      }
    }
  }
  norm_tmp = norm;
  mesh->get_comm()->SumAll(&norm_tmp, &norm, 1);
  std::cout << "Solute mass in matrix: " << 2 * norm << std::endl;

  return err;
}


TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_1)
{
  // (3) diffusion/dispersion in fracture Df
  // (6) Matrix depth, L = 2
  double err = RunTest(1, 0.0, 1e-6, 0.0, 0.0, 2.0);
  CHECK(err < 1.0e-3);
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_2)
{
  // d(a C)/dt + d(u C)/dx - (phi * Dm / (a/2)) (Cm - C) = 0
  double a = 0.001;
  // (2) velocity * aperture = 1e-4 * a
  double u = 1e-4 * a;
  // (4) molecular diffusion coefficient in matrix
  double Dm = 5e-9;
  // (5) normal diffusion coeffcient, kn = phi * Dm / (a/2)
  double kn = 0.2 * Dm / (a / 2);
  double err = RunTest(2, u, 0.0, Dm, kn, 1.0);
  CHECK(err < 1.0e-1);
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_3)
{
  double a = 0.001;
  // (2) velocity * aperture = 1e-4 * a
  double u = 1e-4 * a;
  // (4) molecular diffusion coefficient in matrix
  double Dm = 5e-9;
  // (5) normal diffusion coeffcient, kn = phi * Dm / (a/2)
  double kn = 0.2 * Dm / (a / 2);
  double err = RunTest(3, u, 0.0, Dm, kn, 0.05, 0.5e+5);
  CHECK(err < 0.01);
}

/*
TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_4) {
  double a = 0.02;
  // (2) velocity * aperture = 1e-4 * a
  double u = 1e-4 * a;
  // (4) molecular diffusion coefficient in matrix, Dm = 2e-6
  // (5) normal diffusion coeffcient, kn = phi * Dm / (a/2) = 0.2 * 2e-6 / (a/2)
  double kn = 0.2 * 2e-6 / (a / 2);
  double err = RunTest(4, u, 0.0, 2.0e-6, kn, 2.4e+5, 2.0);
  // CHECK(err < 2.0e-2);
}
*/

/*
TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_5) {
  // d(a C)/dt + d(u C)/dx - a Df d^2(C)/dx^2 - (Dm/2) (Cm - C) = 0
  double err = RunTest(3, 1.0e-5, 1.0e-5, 1.0e-6, 1.0e-6/2);
  CHECK(err < 100);
}
*/
