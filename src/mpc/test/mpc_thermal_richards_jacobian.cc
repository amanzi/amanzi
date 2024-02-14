/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MPC PK
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CompositeVector.hh"
#include "eos_reg.hh"
#include "IO.hh"
#include "FlowEnergy_PK.hh"
#include "MeshFactory.hh"
#include "models_energy_reg.hh"
#include "models_flow_reg.hh"
#include "Operator.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "Richards_PK.hh"
#include "State.hh"
#include "VerboseObject.hh"

void
IndexToSubVector(int n, int nfaces, int ncells, int* i, std::string* kind, int* v)
{
  if (n < nfaces) {
    *i = 0;
    *kind = "face";
    *v = n;
  } else if (n < nfaces + ncells) {
    *i = 0;
    *kind = "cell";
    *v = n - nfaces;
  } else if (n < nfaces + ncells + nfaces) {
    *i = 1;
    *kind = "face";
    *v = n - nfaces - ncells;
  } else {
    *i = 1;
    *kind = "cell";
    *v = n - nfaces - ncells - nfaces;
  }
}


void
ComputeEigenvalues(const Amanzi::WhetStone::DenseMatrix& matrix,
                   Amanzi::WhetStone::DenseVector& real,
                   Amanzi::WhetStone::DenseVector& imag)
{
  Amanzi::WhetStone::DenseMatrix A(matrix);
  int n = A.NumRows();

  int info, ldv(1), lwork = 4 * n;
  double VL, VR, dwork[lwork];

  real.Reshape(n);
  imag.Reshape(n);

  Amanzi::WhetStone::DGEEV_F77("N",
                               "N",
                               &n,
                               A.Values(),
                               &n,
                               real.Values(),
                               imag.Values(),
                               &VL,
                               &ldv,
                               &VR,
                               &ldv,
                               dwork,
                               &lwork,
                               &info);
}


/* ****************************************************************
* Analysis of Jacobian for the FEHM EOS.
* ************************************************************** */
TEST(ENERGY_JACOBIAN)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: Jacobian calculation for thermal Richards" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mpc_thermal_richards_jacobian.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 8, 8);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("flow and energy");
  auto soln = Teuchos::rcp(new TreeVector());
  auto MPC = Teuchos::rcp(new FlowEnergy_PK(pk_tree, plist, S, soln));

  MPC->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  MPC->Initialize();
  S->CheckAllFieldsInitialized();
  MPC->CommitStep(0.0, 1.0, Tags::DEFAULT);

  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("FlowEnergy", *plist));
  WriteStateStatistics(*S, *vo);

  // finite difference Jacobian
  auto u1 = Teuchos::rcp(new TreeVector());
  auto u1f = Teuchos::rcp(new TreeVector());
  auto u1e = Teuchos::rcp(new TreeVector());
  u1->PushBack(u1f);
  u1->PushBack(u1e);
  u1f->SetData(S->GetPtrW<CompositeVector>("pressure", Tags::DEFAULT, ""));
  u1e->SetData(S->GetPtrW<CompositeVector>("temperature", Tags::DEFAULT, ""));

  auto u0 = Teuchos::rcp(new TreeVector(*u1));
  auto f0 = Teuchos::rcp(new TreeVector(*u1));
  auto f1 = Teuchos::rcp(new TreeVector(*u1));

  int ncells, nfaces, nJ, v, i;
  double umax[2], factor, eps(1e-6), t_old(0.0), t_new(1.0);
  std::string kind;

  ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nJ = 2 * (ncells + nfaces);
  WhetStone::DenseMatrix Jfd(nJ, nJ);

  u1->SubVector(0)->NormInf(&umax[0]);
  u1->SubVector(1)->NormInf(&umax[1]);

  Jfd.PutScalar(0.0);
  MPC->FunctionalResidual(t_old, t_new, u0, u1, f0);

  for (int ncol = 0; ncol < nJ; ++ncol) {
    MPC->ChangedSolution();

    IndexToSubVector(ncol, nfaces, ncells, &i, &kind, &v);
    auto& u1_v = *u1->SubVector(i)->Data()->ViewComponent(kind);

    factor = eps * umax[i];
    u1_v[0][v] += factor;
    MPC->FunctionalResidual(t_old, t_new, u0, u1, f1);
    u1_v[0][v] -= factor;

    for (int nrow = 0; nrow < nJ; ++nrow) {
      IndexToSubVector(nrow, nfaces, ncells, &i, &kind, &v);
      auto& f0_v = *f0->SubVector(i)->Data()->ViewComponent(kind);
      auto& f1_v = *f1->SubVector(i)->Data()->ViewComponent(kind);
      Jfd(nrow, ncol) = (f1_v[0][v] - f0_v[0][v]) / factor;
    }
  }

  // numerical Jacobian
  MPC->UpdatePreconditioner(t_old, u1, t_new - t_old);
  auto it = MPC->begin();
  auto pk0 = Teuchos::rcp_dynamic_cast<Flow::Richards_PK>(*it);
  auto A0 = pk0->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->A();

  it++;
  auto pk1 = Teuchos::rcp_dynamic_cast<Energy::Energy_PK>(*it);
  auto A1 = pk1->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->A();

  int num_entries;
  int* indices;
  double* values;

  WhetStone::DenseMatrix Jpk(nJ, nJ);
  Jpk.PutScalar(0.0);

  for (int row = 0; row < nJ / 2; ++row) {
    A0->ExtractMyRowView(row, num_entries, values, indices);
    for (int n = 0; n < num_entries; ++n) {
      int col = indices[n];
      Jpk(row, col) = values[n];
    }
  }

  for (int row = 0; row < nJ / 2; ++row) {
    A1->ExtractMyRowView(row, num_entries, values, indices);
    for (int n = 0; n < num_entries; ++n) {
      int col = indices[n];
      Jpk(nJ / 2 + row, nJ / 2 + col) = values[n];
    }
  }

  // std::cout << Jfd << std::endl;
  // std::cout << Jpk << std::endl;
  auto Jdiff = Jfd - Jpk;
  double jdiff = Jdiff.Norm2();
  double jfd = Jfd.Norm2();
  double jpk = Jpk.Norm2();

  std::cout << "|| Jfd - Jpk || = " << jdiff << ",  || Jfd || = " << jfd << ",  || Jpk || = " << jpk
            << std::endl;
  CHECK(jdiff / jfd < 4e-3);

  WhetStone::DenseVector real, imag;
  WhetStone::DenseMatrix J1(Jpk);
  J1.Inverse();
  auto J2 = J1 * Jfd;
  ComputeEigenvalues(J2, real, imag);

  double tmp, emin(1e+98), emax(0.0);
  for (int n = 0; n < real.NumRows(); ++n) {
    tmp = std::sqrt(real(n) * real(n) + imag(n) * imag(n));
    emin = std::min(emin, tmp);
    emax = std::max(emax, tmp);
  }
  std::cout << "cond(inv(Jpk) * Jfd) = " << emax / emin << std::endl;

  // dFlow / dT in Jacobian changes little in 2-norm but much in spectral norm.
  for (int row = 0; row < nJ / 2; ++row) {
    for (int col = nJ / 2; col < nJ; ++col) Jpk(row, col) = Jfd(row, col);
    for (int col = nJ / 2; col < nJ; ++col) Jpk(col, row) = Jfd(col, row);
  }

  Jdiff = Jfd - Jpk;
  jdiff = Jdiff.Norm2();
  jpk = Jpk.Norm2();

  std::cout << "|| Jfd - Jpk || = " << jdiff << ",  || Jfd || = " << jfd << ",  || Jpk || = " << jpk
            << std::endl;
  CHECK(jdiff / jfd < 4e-3);

  J1 = Jpk;
  J1.Inverse();
  J2 = J1 * Jfd;
  ComputeEigenvalues(J2, real, imag);

  emin = 1e+98;
  emax = 0.0;
  for (int n = 0; n < real.NumRows(); ++n) {
    tmp = std::sqrt(real(n) * real(n) + imag(n) * imag(n));
    emin = std::min(emin, tmp);
    emax = std::max(emax, tmp);
  }
  std::cout << "cond(inv(Jpk) * Jfd) = " << emax / emin << std::endl;
}
