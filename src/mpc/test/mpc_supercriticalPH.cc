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

#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "evaluators_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "models_energy_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"
#include "models_flow_reg.hh"


TEST(MPC_DRIVER_THERMAL_RICHARDS_JACOBIAN)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  auto plist = Teuchos::getParametersFromXmlFile("test/mpc_supercriticalPH.xml");
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");

  auto comm = Amanzi::getDefaultComm();
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 10.0, 3.0, 10, 3);

  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterMesh("domain", mesh);

  auto pk_tree = plist->sublist("cycle driver").sublist("time periods").sublist("TP 0").sublist("PK tree").sublist("transient:mpc1");
  auto soln = Teuchos::rcp(new TreeVector());
  auto mpc = Teuchos::rcp(new FlowEnergyPH_PK(pk_tree, plist, S, soln));

  mpc->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  mpc->Initialize();

  S->CheckAllFieldsInitialized();
  S->InitializeIOFlags();

  // WriteStateStatistics(*S);

  // numerical Jacobian
  double dt0(0.1);
  mpc->UpdatePreconditioner(0.1, soln, dt0);
  auto A = mpc->op_tree_pc()->A();

  int ncells, nfaces;
  ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

  int num_entries, nJ(nfaces + ncells);
  int* indices;
  double* values;

  WhetStone::DenseMatrix Jpk(2 * nJ, 2 * nJ);
  Jpk.PutScalar(0.0);

  for (int row = 0; row < 2 * nJ; ++row) {
    A->ExtractMyRowView(row, num_entries, values, indices);
    for (int n = 0; n < num_entries; ++n) {
      int col = indices[n];
      Jpk(row, col) = values[n];
    }
  }

  // finite difference Jacobian
  auto u1 = soln;

  auto u0 = Teuchos::rcp(new TreeVector(*u1));
  auto f0 = Teuchos::rcp(new TreeVector(*u1));
  auto f1 = Teuchos::rcp(new TreeVector(*u1));
 
  int v;
  double umax, factor, eps(1e-8), t_old(0.0), t_new(dt0);
  std::string kind;

  WhetStone::DenseMatrix Jfd(nJ, nJ), Jsub(nJ, nJ);

  for (int i0 = 0; i0 < 2; ++i0) {
    for (int j0 = 0; j0 < 2; ++j0) {
      u0->SubVector(j0)->NormInf(&umax);
      factor = eps * umax;

      Jfd.PutScalar(0.0);
      mpc->ChangedSolution();
      mpc->FunctionalResidual(t_old, t_new, u0, u1, f0);

      for (int ncol = 0; ncol < nJ; ++ncol) {
        if (ncol < nfaces) {
          kind = "face";
          v = ncol;
        } else {
          kind = "cell";
          v = ncol - nfaces;
        }
        auto& u1_v = *u1->SubVector(j0)->Data()->ViewComponent(kind);

        u1_v[0][v] += factor;
        mpc->ChangedSolution();
        mpc->FunctionalResidual(t_old, t_new, u0, u1, f1);
        u1_v[0][v] -= factor;

        for (int nrow = 0; nrow < nJ; ++nrow) {
          if (nrow < nfaces) {
            kind = "face";
            v = nrow;
          } else {
            kind = "cell";
            v = nrow - nfaces;
          }

          auto& f0_v = *f0->SubVector(i0)->Data()->ViewComponent(kind);
          auto& f1_v = *f1->SubVector(i0)->Data()->ViewComponent(kind);
          Jfd(nrow, ncol) = (f1_v[0][v] - f0_v[0][v]) / factor;
        }
      }

      std::vector<int> idx(nJ), jdx(nJ);
      for (int f = 0; f < nfaces; ++f) {
        idx[f] = 2 * f + i0;
        jdx[f] = 2 * f + j0;
      } 
      for (int c = 0; c < ncells; ++c) {
        idx[nfaces + c] = 2 * nfaces + 2 * c + i0;
        jdx[nfaces + c] = 2 * nfaces + 2 * c + j0;
      }

      for (int i = 0; i < nJ; ++i) {
        for (int j = 0; j < nJ; ++j) { 
          Jsub(i, j) = Jpk(idx[i], jdx[j]);
        }
      }

      // hack: enforce Dirichlet for finite-difference approximation
      if (i0 == j0) {
        for (int i = 0; i < nJ; ++i) {
          if (Jsub(i, i) == 1.0) Jfd(i, i) = 1.0;
        }
      }
  
      auto Jdiff = Jfd - Jsub;
      double norm_diff = Jdiff.Norm2();
      double norm_fd = Jfd.Norm2();
      double norm_pk = Jsub.Norm2();

      std::cout << i0 << j0 << std::scientific 
                << ", || Jfd - Jpk || = " << norm_diff
                << ",  || Jfd || = " << norm_fd 
                << ",  || Jpk || = " << norm_pk << std::endl;
      CHECK(norm_diff < 1e-3 * norm_pk);
      int s1(0), s2(266);
      // if (i0 == 1 && j0 == 0) std::cout << Jfd.SubMatrix(s1, s1 + 5, s2, s2 + 8) << std::endl;
      // if (i0 == 1 && j0 == 0) std::cout << Jsub.SubMatrix(s1, s1 + 5, s2, s2 + 8) << std::endl;
      if (i0 == 1 && j0 == 0) {
        // std::cout << Jfd << std::endl;
        // std::cout << Jsub << std::endl;
      }
    }
  }
}


TEST(MPC_DRIVER_THERMAL_RICHARDS_HIGH_PH)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  auto plist = Teuchos::getParametersFromXmlFile("test/mpc_supercriticalPH.xml");
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");

  auto comm = Amanzi::getDefaultComm();
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 100.0, 4.0, 200, 2);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();
}

