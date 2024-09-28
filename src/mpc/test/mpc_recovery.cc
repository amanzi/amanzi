/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MPC PK

  Tests for failure/recovery
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
#include "bilinear_form_reg.hh"
#include "CompositeVector.hh"
#include "eos_reg.hh"
#include "evaluators_flow_reg.hh"
#include "evaluators_multiphase_reg.hh"
#include "IO.hh"
#include "MeshFactory.hh"
#include "MeshExtractedManifold.hh"
#include "models_energy_reg.hh"
#include "models_flow_reg.hh"
#include "models_multiphase_reg.hh"
#include "Operator.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mechanics_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_multiphase_reg.hh"
#include "Richards_PK.hh"
#include "State.hh"
#include "VerboseObject.hh"

// Amanzi::MPC
#include "mpc_utils.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

/* ****************************************************************
* Analysis of Failure/Recovery 
* ************************************************************** */
template <class PK>
void
Run(const std::string& xmlFileName, int dim, const std::vector<double>& dt, int iTP, int icase = -1)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTEST: MPC failure and recovery" << std::endl;

  // read parameter list
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(dim, region_list, *comm));

  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  Teuchos::ParameterList& aux = mesh_list->sublist("unstructured");
  bool has_submesh = aux.isSublist("submesh");
  if (has_submesh) {
    mesh_list->set<bool>("request edges", true);
    mesh_list->set<bool>("request faces", true);
  }

  MeshFactory meshfactory(comm, gm, mesh_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  Teuchos::RCP<const Mesh> mesh, submesh;
  if (aux.isSublist("generate mesh")) {
    Teuchos::ParameterList tmp = aux.sublist("generate mesh");
    mesh = meshfactory.create(tmp);
  } else {
    std::string filename = aux.sublist("read mesh file").get<std::string>("file");
    mesh = meshfactory.create(filename);
  }

  // create additional mesh for fracture
  if (has_submesh) {
    auto names = aux.sublist("submesh").get<Teuchos::Array<std::string>>("regions").toVector();
    auto submesh_fr = Teuchos::rcp(
      new MeshExtractedManifold(mesh, names[0], AmanziMesh::Entity_kind::FACE, comm, gm, plist));
    submesh =
      Teuchos::rcp(new Mesh(submesh_fr, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), plist));
  }

  // create a simple state and populate it
  std::vector<Teuchos::RCP<State>> S(2);
  std::vector<Teuchos::RCP<PK>> pk(2);
  std::vector<Teuchos::RCP<TreeVector>> soln(2);

  for (int i = 0; i < 2; ++i) {
    std::string sTP = "TP " + std::to_string(iTP);
    Teuchos::ParameterList tmp =
      plist->sublist("cycle driver").sublist("time periods").sublist(sTP).sublist("PK tree");
    std::string name = tmp.begin()->first;
    Teuchos::ParameterList pk_tree = tmp.sublist(name);

    if (i == 1)
      plist->sublist("PKs")
        .sublist(name)
        .sublist("time integrator")
        .sublist("BDF1")
        .set<int>("report failure on step", 1);

    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
    state_list.sublist("verbose object").set<std::string>("verbosity level", "none");
    S[i] = Teuchos::rcp(new State(state_list));

    // special cases
    if (icase == 5) {
      S[i]->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(submesh));

      Key key("mass_density_gas");
      S[i]
        ->Require<CompositeVector, CompositeVectorSpace>(key, Tags::DEFAULT, key)
        .SetMesh(submesh)
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    } else if (icase == 6) {
      S[i]->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

      Key key("mass_density_gas");
      S[i]
        ->Require<CompositeVector, CompositeVectorSpace>(key, Tags::DEFAULT, key)
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    } else {
      S[i]->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));
      if (has_submesh) S[i]->RegisterMesh("fracture", Teuchos::rcp_const_cast<Mesh>(submesh));
    }

    soln[i] = Teuchos::rcp(new TreeVector());
    pk[i] = Teuchos::rcp(new PK(pk_tree, plist, S[i], soln[i]));

    pk[i]->Setup();
    S[i]->Setup();
    S[i]->InitializeFields();
    S[i]->InitializeEvaluators();

    pk[i]->Initialize();
    S[i]->CheckAllFieldsInitialized();
    pk[i]->CommitStep(0.0, 1.0, Tags::DEFAULT);
    pk[i]->CalculateDiagnostics(Tags::DEFAULT);

    // WriteStateStatistics(*S[i]);
  }

  // perform given steps
  int ndts = dt.size();
  double told(0.0), tnew;
  for (int n = 0; n < ndts; ++n) {
    tnew = told + dt[n];
    pk[0]->AdvanceStep(told, tnew, true);
    pk[0]->CommitStep(told, tnew, Tags::DEFAULT);
    told = tnew;
  }

  // perform failure on the second step
  std::cout << std::endl << std::endl;
  std::vector<double> dt_new({ dt[0], dt[1], dt[1], dt[2] });

  for (int n = 0; n < ndts + 1; ++n) {
    tnew = told + dt_new[n];
    bool fail = pk[1]->AdvanceStep(told, tnew, true);
    if (!fail) {
      pk[1]->CommitStep(told, tnew, Tags::DEFAULT);
      told = tnew;
    }
  }

  // compare states
  for (int i = 0; i < 2; ++i) { pk[i]->CalculateDiagnostics(Tags::DEFAULT); }

  double err, fnorm;
  for (auto r = S[0]->data_begin(); r != S[0]->data_end(); ++r) {
    std::string name = r->first;

    for (const auto& r : S[0]->GetRecordSet(name)) {
      if (r.second->ValidType<CompositeVector>()) {
        const auto& f0 = r.second->Get<CompositeVector>();
        const auto& f1 = S[1]->Get<CompositeVector>(name, Tags::DEFAULT);

        auto df(f1);
        df.Norm2(&fnorm);
        df.Update(1.0, f0, -1.0);
        df.Norm2(&err);
        if (fnorm != 0.0) err /= fnorm;

        printf("%33s:  err=%12.6g   norm=%12.6g\n", name.c_str(), err, fnorm);
        CHECK(err < 1e-14);
      }
    }
  }

  // WriteStateStatistics(*S[0]);
  // WriteStateStatistics(*S[1]);
}

TEST(MPC_RECOVERY_COUPLED_THERMAL_FLOW)
{
  CreateApertureFile(144, 300.0);

  ::Run<FlowEnergyMatrixFracture_PK>(
    "test/mpc_coupled_thermal_flow_richards.xml", 3, { 10.0, 10.0, 10.0 }, 1);
    // round-off errors are observed between two states/two runs
    // "test/mpc_coupled_thermal_flow_richards.xml", 3, { 0.1, 0.1, 0.1 }, 1);
}

TEST(MPC_RECOVERY_FLOW_RICHARDS)
{
  ::Run<Flow::Richards_PK>("test/mpc_flow.xml", 2, { 15.0, 15.0, 15.0 }, 0);
}

TEST(MPC_RECOVERY_THERMAL_FLOW)
{
  ::Run<FlowEnergy_PK>("test/mpc_thermal_richards.xml", 2, { 1.0e+4, 1.0e+4, 1.0e+4 }, 0);
}

TEST(MPC_RECOVERY_COUPLED_FLOW)
{
  ::Run<FlowMatrixFracture_PK>("test/mpc_coupled_flow.xml", 3, { 1.0e+6, 1.0e+6, 1.0e+6 }, 0);
}

TEST(MPC_RECOVERY_MULTIPHASE_FRACTURES)
{
  ::Run<Multiphase::Multiphase_PK>(
    "test/mpc_multiphase_fractures.xml", 3, { 2.0e+11, 2.0e+11, 2.0e+11 }, 0, 5);
}

TEST(MPC_RECOVERY_MULTIPHASE_THERMAL_3D)
{
  ::Run<Multiphase::Multiphase_PK>(
    "test/mpc_multiphase_thermal.xml", 3, { 2.0e+7, 2.0e+7, 2.0e+7 }, 0, 6);
}
