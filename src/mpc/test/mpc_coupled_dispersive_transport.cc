/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

  Solution of PDE:

    d(a C)/dt + d(u C)/dx - Df d^2(C)/dx^2 \sum_i (phi * Dm_i) d(Cm_i)/dy = 0
 
  where a is aperture, u is velocity [m2/s], phi is matrix porosity,
  Dm is molecular diffusion, Df is dispersion/diffusion in fracture and Cm_i is 
  matrix concentration on the i-th side of the fracture. 

  Case 1: Dm=0 and u=0
  Case 2: Df=0 and constant BC on the fracture inlet
  Case 3: Df=0 and exponential decaying BC on the fracture inlet
  // Case 4: Df=0 and pulse BC on the fracture inlet
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
#include "eos_reg.hh"
#include "evaluators_mpc_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "models_transport_reg.hh"
#include "NumericalIntegration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"

// MPC
#include "mpc_utils.hh"

double
Case3(double a, double b, double c, double t)
{
  return a * pow(t, -1.5) * std::exp(-b / t + c * t);
}


double
RunTest(int icase, double u_f, double mol_diff_f, double mol_diff_m, double L, double tend = 1e+5)
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

  // coefficient for solute diffusion to matrix
  plist->sublist("state")
    .sublist("evaluators")
    .sublist("fracture-solute_diffusion_to_matrix")
    .set<double>("molecular diffusion", mol_diff_m);

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
      .set<std::string>("formula", "if(t > 0.5, 0.0, 1.0)");

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
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 3.0 - L, 7.0, 1.0, 3.0 + L, 56, 1, 32);

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

  const auto& tcc_m =
    *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  double cmin(1e+99), cmax(-1e+99);
  for (int c = 0; c < tcc_m.MyLength(); ++c) {
    const auto& xc = mesh->getCellCentroid(c);
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
  double phi = (*S->Get<CompositeVector>("porosity").ViewComponent("cell"))[0][0];

  double t(tend), Df(mol_diff_f), Dm(mol_diff_m);

  for (int c = 0; c < tcc_f.MyLength(); ++c) {
    const auto& xc = mesh_fracture->getCellCentroid(c);
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
        double tmp2 = std::erfc(xc[0] * phi * std::sqrt(Dm / (u_f * b * Tm)));
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
  if (icase != 4) {
    mesh->get_comm()->SumAll(&err_tmp, &err, 1);
    mesh->get_comm()->SumAll(&norm_tmp, &norm, 1);
    err /= tcc_f.GlobalLength();
    norm /= tcc_f.GlobalLength();
    std::cout << "Mean error in fracture: " << err << " solution norm=" << norm << std::endl;

  } else {
    auto data = getObservations("observation.out", 9);

    double t0, t, a, tmp0, tmp1, tmp2, norm(0.0);
    t0 = 3 * b / u_f;
    a = phi * sqrt(Dm) / b;
    tmp0 = a * sqrt(1.0 / t0 / M_PI);

    for (auto d : data) {
      t = d.first;
      if (t > t0) {
        tmp1 = t / t0 - 1.0;
        tmp2 = tmp0 * std::pow(tmp1, -1.5) * std::exp(-a * a * t0 / tmp1);
        err += std::pow(d.second - tmp2, 2.0);
        norm += std::pow(tmp2, 2.0);
      }
    }
    err = std::sqrt(err) / data.size();
    norm = std::sqrt(err) / data.size();
    std::cout << "Mean error in observations: " << err << " norm=" << norm << std::endl;
  }

  if (fmin < 1.0e+98) CHECK_CLOSE(fmin, fmax, 1e-12);


  // mass of solute in matrix
  // err = 0.0;
  norm = 0.0;
  for (int c = 0; c < tcc_m.MyLength(); ++c) {
    const auto& xc = mesh->getCellCentroid(c);
    double vol = mesh->getCellVolume(c);

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
  mesh->getComm()->SumAll(&norm_tmp, &norm, 1);
  std::cout << "Solute mass in matrix: " << 2 * norm << std::endl;

  return err;
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_1)
{
  double L = 2.0; // matrix depth
  double Df = 1e-6;
  double Dm = 0.0;
  double err = RunTest(1, 0.0, Df, Dm, L);
  CHECK(err < 1.0e-3);
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_2)
{
  double a = 0.001;
  double u = 1e-4 * a;
  double Dm = 5e-9;
  double err = RunTest(2, u, 0.0, Dm, 0.25);
  CHECK(err < 1.0e-1);
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_3)
{
  double a = 0.001;
  double u = 1e-4 * a;
  double Dm = 5e-9;
  double err = RunTest(3, u, 0.0, Dm, 0.05, 0.5e+5);
  CHECK(err < 0.01);
}

TEST(MPC_DIFFUSIVE_TRANSPORT_MATRIX_FRACTURE_4)
{
  double a = 0.001;
  double u = 1e-4 * a;
  double Dm = 5e-10;
  double err = RunTest(4, u, 0.0, Dm, 0.05, 1e+5);
  CHECK(err < 2.0e-2);
}
