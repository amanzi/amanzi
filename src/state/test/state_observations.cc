/*

  Copyright 2010-201x held jointly by LANS/LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "CompositeVector.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "Observable.hh"
#include "UnstructuredObservations.hh"

using namespace Amanzi;

//
//  shamelessly stole from SO
//
// https://stackoverflow.com/questions/6163611/compare-two-files
//
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>

bool compareFiles(const std::string& p1, const std::string& p2) {
  std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; //file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false; //size mismatch
  }

  //seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}


struct obs_test {
public:
  obs_test() {
    auto comm = Amanzi::getDefaultComm();

    // create geometric model
    Teuchos::ParameterList region_list("regions");

    auto& top = region_list.sublist("top face").sublist("region: plane");
    top.set<Teuchos::Array<double>>("point", std::vector<double>{1,1,1});
    top.set<Teuchos::Array<double>>("normal", std::vector<double>{0,0,1});

    region_list.sublist("all").sublist("region: all");

    auto& middle = region_list.sublist("middle").sublist("region: box");
    middle.set<Teuchos::Array<double>>("low coordinate", std::vector<double>{-0.3,-0.3,-0.3});
    middle.set<Teuchos::Array<double>>("high coordinate", std::vector<double>{0.3,0.3,0.3});

    auto& one_side = region_list.sublist("one side volume").sublist("region: box");
    one_side.set<Teuchos::Array<double>>("low coordinate", std::vector<double>{-1, -1, -1});
    one_side.set<Teuchos::Array<double>>("high coordinate", std::vector<double>{-0.3, 1, 1});

    auto& one_side_side = region_list.sublist("one side side").sublist("region: box");
    one_side_side.set<Teuchos::Array<double>>("low coordinate", std::vector<double>{-.35, -1, -1});
    one_side_side.set<Teuchos::Array<double>>("high coordinate", std::vector<double>{-0.3, 1, 1});

    Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

    AmanziMesh::MeshFactory meshfactory(comm,gm);
    Teuchos::RCP<AmanziMesh::Mesh> mesh = meshfactory.create(-1,-1,-1,1,1,1,3,3,3);

    Teuchos::ParameterList state_list("state");

    S = Teuchos::rcp(new State(state_list));
    S->RegisterMesh("domain", mesh);
    S->RequireField("constant")->SetMesh(mesh)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireField("linear")->SetMesh(mesh)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireField("id")->SetMesh(mesh)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireField("flux")->SetMesh(mesh)->SetGhosted(false)
      ->SetComponent("face", AmanziMesh::FACE, 1);
    S->set_time(0.);
    S->set_cycle(0);

    S->Setup();
    S->GetFieldData("constant", "state")->PutScalar(2.0);
    S->GetFieldData("linear", "state")->PutScalar(0.);


    Epetra_MultiVector& flux_f = *S->GetFieldData("flux", "state")
      ->ViewComponent("face", false);
    AmanziGeometry::Point plus_xz(1.0, 0.0, 1.0);
    for (int f=0; f!=flux_f.MyLength(); ++f) {
      flux_f[0][f] = mesh->face_normal(f) * plus_xz;
    }

    Epetra_MultiVector& id_c = *S->GetFieldData("id", "state")
      ->ViewComponent("cell", false);
    auto& cell_map = S->GetMesh("domain")->map(AmanziMesh::CELL, false);
    for (int c=0; c!=id_c.MyLength(); ++c) {
      id_c[0][c] = cell_map.GID(c);
    }
  }

  void advance(double dt) {
    S->advance_time(dt);
    S->GetFieldData("linear", "state")->PutScalar(S->time() * 0.1);
    S->advance_cycle();
  }

public:
  Teuchos::RCP<State> S;
};

SUITE(STATE_OBSERVATIONS) {

TEST_FIXTURE(obs_test, Assumptions) {
  int num_cells = S->GetMesh("domain")
    ->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int num_cells_total = 0;
  S->GetMesh("domain")->get_comm()->SumAll(&num_cells, &num_cells_total, 1);
  CHECK_EQUAL(27, num_cells_total);

  int faces_middle = S->GetMesh("domain")
    ->get_set_size("middle", AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(0, faces_middle);

  int cells_all = S->GetMesh("domain")
    ->get_set_size("all", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int cells_all_total = 0;
  S->GetMesh("domain")->get_comm()->SumAll(&cells_all, &cells_all_total, 1);
  CHECK_EQUAL(27, cells_all_total);
}

TEST_FIXTURE(obs_test, ObservePoint) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "constant");
  obs_list.set<std::string>("region", "middle");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "point");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(2.0, observation[0], 1.e-10);
}

TEST_FIXTURE(obs_test, ObserveIntensiveIntegral) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "constant");
  obs_list.set<std::string>("region", "all");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "integral");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(2.0* (2*2*2), observation[0], 1.e-10); // total volume
}


TEST_FIXTURE(obs_test, ObserveExtensiveIntegral) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "constant");
  obs_list.set<std::string>("region", "all");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "extensive integral");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(2.0*27, observation[0], 1.e-10); // number of cells
}


TEST_FIXTURE(obs_test, ObserveAverage) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "constant");
  obs_list.set<std::string>("region", "all");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "average");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(2.0, observation[0], 1.e-10); // vol weighted average
}

TEST_FIXTURE(obs_test, ObserveMin) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "id");
  obs_list.set<std::string>("region", "all");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "minimum");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(0.0, observation[0], 1.e-10); // smallest id is 0
}

TEST_FIXTURE(obs_test, ObserveMax) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "id");
  obs_list.set<std::string>("region", "all");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "maximum");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);

  CHECK_CLOSE(26, observation[0], 1.e-10); // biggest id is num cells - 1
}


TEST_FIXTURE(obs_test, Face) {
  // integrate an observable
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "flux");
  obs_list.set<std::string>("region", "top face");
  obs_list.set<std::string>("location name", "face");
  obs_list.set<std::string>("functional", "extensive integral");
  obs_list.set("direction normalized flux", true);

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(4.0, observation[0], 1.e-10); // flux is given by face area
}


TEST_FIXTURE(obs_test, Face_NORMALIZED_REL_VOLUME) {
  // direction nomralized flux relative to region allows normalizing in an
  // outward-normal relative to a volumetric region.
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("variable", "flux");
  obs_list.set<std::string>("region", "one side side");
  obs_list.set<std::string>("location name", "face");
  obs_list.set<std::string>("functional", "extensive integral");
  obs_list.set("direction normalized flux", true);
  obs_list.set("direction normalized flux relative to region", "one side volume");

  Observable obs(obs_list);
  obs.Setup(S.ptr());
  obs.FinalizeStructure(S.ptr());
  CHECK_EQUAL(1, obs.get_num_vectors());

  std::vector<double> observation(1, Observable::nan);
  obs.Update(S.ptr(), observation, 0);
  CHECK_CLOSE(4.0, observation[0], 1.e-10); // flux is given by face area
}



TEST_FIXTURE(obs_test, FileOne) {
  //  one observation in a file
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("observation output filename", "obs1.dat");
  obs_list.set<Teuchos::Array<int>>("cycles", std::vector<int>{0,1});
  obs_list.set<std::string>("variable", "linear");
  obs_list.set<std::string>("region", "all");
  obs_list.set<std::string>("location name", "cell");
  obs_list.set<std::string>("functional", "average");

  {
    UnstructuredObservations obs(obs_list);
    obs.Setup(S.ptr());
    obs.MakeObservations(S.ptr());
    advance(1.0);
    obs.MakeObservations(S.ptr());
  }

  // times: 0, 1
  // values: 0, .1
  CHECK(compareFiles("obs1.dat", "test/obs1.dat.gold"));
}


TEST_FIXTURE(obs_test, FileTwo) {
  //  one observation in a file
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("observation output filename", "obs2.dat");
  obs_list.set<Teuchos::Array<int>>("cycles", std::vector<int>{0,1});

  auto& obsA_list = obs_list.sublist("observed quantities").sublist("obsA");
  obsA_list.set<std::string>("variable", "linear");
  obsA_list.set<std::string>("region", "all");
  obsA_list.set<std::string>("location name", "cell");
  obsA_list.set<std::string>("functional", "average");

  auto& obsB_list = obs_list.sublist("observed quantities").sublist("obsB");
  obsB_list.set<std::string>("variable", "constant");
  obsB_list.set<std::string>("region", "all");
  obsB_list.set<std::string>("location name", "cell");
  obsB_list.set<std::string>("functional", "average");

  {
    UnstructuredObservations obs(obs_list);
    obs.Setup(S.ptr());
    obs.MakeObservations(S.ptr());
    advance(1.0);
    obs.MakeObservations(S.ptr());
  }

  // times: 0, 1
  // valuesA: 0, .1
  // valuesB: 2, 2
  CHECK(compareFiles("obs2.dat", "test/obs2.dat.gold"));
}



TEST_FIXTURE(obs_test, TimeIntegrated) {
  //  one observation in a file
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("observation output filename", "obs3.dat");
  obs_list.set<Teuchos::Array<double>>("times", std::vector<double>{0,1});

  auto& obsA_list = obs_list.sublist("observed quantities").sublist("obsA");
  obsA_list.set<std::string>("variable", "linear");
  obsA_list.set<std::string>("region", "all");
  obsA_list.set<std::string>("location name", "cell");
  obsA_list.set<std::string>("functional", "average");
  obsA_list.set("time integrated", true);

  auto& obsB_list = obs_list.sublist("observed quantities").sublist("obsB");
  obsB_list.set<std::string>("variable", "constant");
  obsB_list.set<std::string>("region", "all");
  obsB_list.set<std::string>("location name", "cell");
  obsB_list.set<std::string>("functional", "average");

  {
    UnstructuredObservations obs(obs_list);
    obs.Setup(S.ptr());
    obs.MakeObservations(S.ptr());
    advance(0.5);
    obs.MakeObservations(S.ptr());
    advance(0.5);
    obs.MakeObservations(S.ptr());
  }

  // times: 0, 1
  // valuesA: 0, 0.5*0.05 + 0.5*0.1  (0.075)
  // valuesB: 2, 2
  CHECK(compareFiles("obs3.dat", "test/obs3.dat.gold"));
}

TEST_FIXTURE(obs_test, WritesNaN) {
  // integrate an observable
  //  one observation in a file
  Teuchos::ParameterList obs_list("my obs");
  obs_list.set<std::string>("observation output filename", "obs4.dat");
  obs_list.set<Teuchos::Array<double>>("times", std::vector<double>{0,1});

  auto& obsA_list = obs_list.sublist("observed quantities").sublist("obsA");
  obsA_list.set<std::string>("variable", "flux");
  obsA_list.set<std::string>("region", "middle");
  obsA_list.set<std::string>("location name", "face"); // no faces in middle region!
  obsA_list.set<std::string>("functional", "extensive integral");

  {
    UnstructuredObservations obs(obs_list);
    obs.Setup(S.ptr());
    obs.MakeObservations(S.ptr());
    advance(0.5);
    obs.MakeObservations(S.ptr());
    advance(0.5);
    obs.MakeObservations(S.ptr());
  }

  // times: 0, 1
  // valuesA: NaN NaN
  CHECK(compareFiles("obs4.dat", "test/obs4.dat.gold"));
}


}
