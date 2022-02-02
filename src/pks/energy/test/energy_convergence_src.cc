/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "AmanziComm.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "EvaluatorPrimary.hh"
#include "GMVMesh.hh"
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Energy
#include "Analytic00.hh"
#include "EnergyOnePhase_PK.hh"
#include "EnthalpyEvaluator.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Energy;

namespace Amanzi {

class TestEnthalpyEvaluator : public EnthalpyEvaluator {
 public:
  explicit TestEnthalpyEvaluator(Teuchos::ParameterList& plist) :
      EnthalpyEvaluator(plist) {};

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new TestEnthalpyEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Evaluate_(
          const State& S, const std::vector<CompositeVector*>& results) override {
    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const auto& temp_c = *S.Get<CompositeVector>("temperature").ViewComponent(*comp);
      auto& result_c = *results[0]->ViewComponent(*comp);
      int ncomp = results[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_c[0][i] = temp_c[0][i];
      }
    }
  }

  virtual void EvaluatePartialDerivative_(
          const State& S, const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& results) override {
    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      auto& result_c = *results[0]->ViewComponent(*comp);
      int ncomp = results[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_c[0][i] = 1.0;
      }
    }
  }

  using EnthalpyEvaluator::EvaluateFieldSingle;
  double EvaluateFieldSingle(const Teuchos::Ptr<State>& S, int c, double T) { return 0.0; }
};

}  // namespace Amanzi


TEST(ENERGY_CONVERGENCE_SRC) {
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout <<"Convergence analysis on three random meshes" << std::endl;

  std::string xmlFileName = "test/energy_convergence_src.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  int nmeshes = plist->get<int>("number of meshes", 1);
  std::vector<double> h, error;

  Teuchos::RCP<VerboseObject> vo_ = Teuchos::rcp(new VerboseObject("", *plist));

  int nx(20);
  double dt(5.0);
  for (int n = 0; n < nmeshes; n++, nx *= 2) {
    dt /= 2.0;

    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));
    
    Preference pref;
    pref.clear();
    pref.push_back(Framework::MSTK);

    MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(pref);
    Teuchos::RCP<const Mesh> mesh;
    if (n == 0) {
      mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 20, 5);
      // mesh = meshfactory.create("test/random_mesh1.exo");
    } else if (n == 1) {
      mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 40, 5);
      // mesh = meshfactory.create("test/random_mesh2.exo");
    } else if (n == 2) {
      mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 80, 5);
      // mesh = meshfactory.create("test/random_mesh3.exo");
    }

    // create a simple state and populate it
    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::ParameterList pk_tree = plist->sublist("PKs").sublist("energy");
    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Teuchos::RCP<EnergyOnePhase_PK> EPK = Teuchos::rcp(new EnergyOnePhase_PK(pk_tree, plist, S, soln));

    // overwrite enthalpy with a different model
    Teuchos::ParameterList ev_list;
    ev_list.set<std::string>("enthalpy key", "enthalpy")
           .set<bool>("include work term", false)
           .set<std::string>("tag", "");
    ev_list.setName("enthalpy");

    S->Require<CompositeVector, CompositeVectorSpace>("enthalpy", Tags::DEFAULT, "enthalpy")
      .SetMesh(mesh)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);

    S->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        "enthalpy", Tags::DEFAULT, "temperature", Tags::DEFAULT, "enthalpy");

    auto enthalpy = Teuchos::rcp(new TestEnthalpyEvaluator(ev_list));
    S->SetEvaluator("enthalpy", enthalpy);

    EPK->Setup();
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    EPK->Initialize();
    S->CheckAllFieldsInitialized();
       
    // constant time stepping 
    int itrs(0);
    double t(0.0), t1(100), dt_next;
    while (t < t1) {
      // swap conserved quntity (no backup, we check dt_next instead)
      const auto& e = S->Get<CompositeVector>("energy");
      auto& e_prev = S->GetW<CompositeVector>("prev_energy", Tags::DEFAULT, "thermal");
      e_prev = e;

      if (itrs == 0) {
        Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
        udot->PutScalar(0.0);
        EPK->bdf1_dae()->SetInitialState(t, soln, udot);
        EPK->UpdatePreconditioner(t, soln, dt);
      }

      EPK->bdf1_dae()->TimeStep(dt, dt_next, soln);
      CHECK(dt_next >= dt);
      EPK->bdf1_dae()->CommitSolution(dt, soln);
      Teuchos::rcp_static_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace> >(
          S->GetEvaluatorPtr("temperature"))->SetChanged();

      t += dt;
      itrs++;
    }

    EPK->CommitStep(0.0, 1.0, Tags::DEFAULT);

    
    // calculate errors
    auto temp = S->GetPtr<CompositeVector>("temperature");
    Analytic00 ana(temp, mesh);

    double l2_norm, l2_err, inf_err;  // error checks
    ana.ComputeCellError(*temp->ViewComponent("cell"), t1, l2_norm, l2_err, inf_err);

    h.push_back(1.0 / nx);
    error.push_back(l2_err);

    printf("mesh=%d bdf1_steps=%3d  L2_temp_err=%7.3e L2_temp=%7.3e\n", n, itrs, l2_err, l2_norm);
    CHECK(l2_err < 8e-1);
    //WriteStateStatistics(*S, *vo_);
    
    // save solution
    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(*temp->ViewComponent("cell"), 0, "temperature");
    GMV::close_data_file();
  }

  
  // check convergence rate
  double l2_rate = Amanzi::Utils::bestLSfit(h, error);
  printf("convergence rate: %10.2f\n", l2_rate);
  CHECK(l2_rate > 0.91);
}


