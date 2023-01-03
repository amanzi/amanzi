/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "CompositeVector.hh"
#include "IO.hh"
#include "MeshFrameworkColumn.hh"
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

bool
compareFiles(const std::string& p1, const std::string& p2)
{
  std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

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
  using CV = CompositeVector;
  using CVS = CompositeVectorSpace;

 public:
  obs_test()
  {
    auto comm = Amanzi::getDefaultComm();

    // create geometric model
    Teuchos::ParameterList region_list("regions");

    auto& top = region_list.sublist("top face").sublist("region: plane");
    top.set<Teuchos::Array<double>>("point", std::vector<double>{ 1, 1, 1 });
    top.set<Teuchos::Array<double>>("normal", std::vector<double>{ 0, 0, 1 });

    region_list.sublist("all").sublist("region: all");

    auto& middle = region_list.sublist("middle").sublist("region: box");
    middle.set<Teuchos::Array<double>>("low coordinate", std::vector<double>{ -0.3, -0.3, -0.3 });
    middle.set<Teuchos::Array<double>>("high coordinate", std::vector<double>{ 0.3, 0.3, 0.3 });

    auto& one_side = region_list.sublist("one side volume").sublist("region: box");
    one_side.set<Teuchos::Array<double>>("low coordinate", std::vector<double>{ -1, -1, -1 });
    one_side.set<Teuchos::Array<double>>("high coordinate", std::vector<double>{ -0.3, 1, 1 });

    auto& one_side_side = region_list.sublist("one side side").sublist("region: box");
    one_side_side.set<Teuchos::Array<double>>("low coordinate",
                                              std::vector<double>{ -.35, -1, -1 });
    one_side_side.set<Teuchos::Array<double>>("high coordinate", std::vector<double>{ -0.3, 1, 1 });

    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

    auto plist = Teuchos::rcp(new Teuchos::ParameterList("mesh factory"));
    plist->sublist("unstructured").sublist("expert").set<std::string>("partitioner", "zoltan_rcb");
    AmanziMesh::MeshFactory meshfactory(comm, gm, plist);
    Teuchos::RCP<AmanziMesh::Mesh> mesh = meshfactory.create(-1, -1, -1, 1, 1, 1, 3, 3, 3);

    Teuchos::ParameterList state_list("state");
    state_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");

    S = Teuchos::rcp(new State(state_list));
    S->RegisterMesh("domain", mesh);

    S->Require<CompositeVector, CompositeVectorSpace>("constant", Tags::DEFAULT, "my_password")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->Require<CompositeVector, CompositeVectorSpace>("linear", Tags::DEFAULT, "my_password")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->Require<CompositeVector, CompositeVectorSpace>("id", Tags::DEFAULT, "my_password")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->Require<CompositeVector, CompositeVectorSpace>("flux", Tags::DEFAULT, "my_password")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    S->Require<CompositeVector, CompositeVectorSpace>("multi_dof", Tags::DEFAULT, "my_password")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 3);
    S->Setup();

    S->set_time(0.0);
    S->set_cycle(0);
  }

  void setup()
  {
    S->Setup();
    S->GetW<CompositeVector>("constant", Tags::DEFAULT, "my_password").PutScalar(2.0);
    S->GetW<CompositeVector>("linear", Tags::DEFAULT, "my_password").PutScalar(0.0);

    (*S->GetW<CV>("multi_dof", Tags::DEFAULT, "my_password").ViewComponent("cell"))(0)->PutScalar(
      0.0);
    (*S->GetW<CV>("multi_dof", Tags::DEFAULT, "my_password").ViewComponent("cell"))(1)->PutScalar(
      1.0);
    (*S->GetW<CV>("multi_dof", Tags::DEFAULT, "my_password").ViewComponent("cell"))(2)->PutScalar(
      2.0);

    auto mesh = S->GetMesh("domain");
    Epetra_MultiVector& flux_f =
      *S->GetW<CV>("flux", Tags::DEFAULT, "my_password").ViewComponent("face");
    AmanziGeometry::Point plus_xz(1.0, 0.0, 1.0);

    for (int f = 0; f != flux_f.MyLength(); ++f) { flux_f[0][f] = mesh->getFaceNormal(f) * plus_xz; }

    Epetra_MultiVector& id_c =
      *S->GetW<CV>("id", Tags::DEFAULT, "my_password").ViewComponent("cell");
    auto& cell_map = S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false);

    for (int c = 0; c != id_c.MyLength(); ++c) { id_c[0][c] = cell_map.GID(c); }
  }

  void advance(double dt)
  {
    S->advance_time(dt);
    S->GetW<CV>("linear", Tags::DEFAULT, "my_password").PutScalar(S->get_time() * 0.1);
    S->advance_cycle();
  }

 public:
  Teuchos::RCP<State> S;
};


struct obs_domain_set_test : public obs_test {
  obs_domain_set_test() : obs_test()
  {
    // create the surface mesh
    auto parent = S->GetMesh("domain");

    auto plist = Teuchos::rcp(new Teuchos::ParameterList("mesh factory"));
    plist->sublist("unstructured").sublist("expert").set<std::string>("partitioner", "zoltan_rcb");
    plist->set<bool>("request faces",true); 
    plist->set<bool>("request edges",false); 
    AmanziMesh::MeshFactory fac(parent->getComm(), parent->getGeometricModel(), plist);
    auto surface_mesh = fac.create(parent, { "top face" }, AmanziMesh::Entity_kind::FACE, true);
    S->RegisterMesh("surface", surface_mesh);

    // create domain set
    assert(false); 
    //parent->buildColumns();
    std::vector<std::string> cols;
    assert(false); 
    //for (int i = 0; i != parent->columns(); ++i) {
    //  cols.emplace_back(std::to_string(surface_mesh->getMap(AmanziMesh::Entity_kind::CELL,false).GID(i)));
    //}
    auto domain_set =
      Teuchos::rcp(new AmanziMesh::DomainSet("column", S->GetMesh("surface"), cols));
    S->RegisterDomainSet("column", domain_set);

    // create subdomain meshes
    int i = 0;
    for (auto& ds : *domain_set) {
      auto parent_list = Teuchos::rcp(new Teuchos::ParameterList(*parent->getParameterList()));
      assert(false); 
      //Teuchos::RCP<AmanziMesh::Mesh> col_mesh =
      //  AmanziMesh::createColumnMesh(parent, i, parent_list);
      //S->RegisterMesh(ds, col_mesh);
      i++;
    }

    for (auto& ds : *domain_set) {
      S->Require<CompositeVector, CompositeVectorSpace>(
         Keys::getKey(ds, "variable"), Tags::DEFAULT, "my_password")
        .SetMesh(S->GetMesh(ds))
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  void setup_domain_set()
  {
    obs_test::setup();
    auto ds = S->GetDomainSet("column");
    for (auto& dname : *ds) {
      int index = Keys::getDomainSetIndex<int>(dname);
      S->GetW<CompositeVector>(Keys::getKey(dname, "variable"), Tags::DEFAULT, "my_password")
        .PutScalar(index);
    }
  }
};


SUITE(STATE_OBSERVATIONS)
{
  TEST_FIXTURE(obs_test, Assumptions)
  {
    setup();
    int num_cells =
      S->GetMesh("domain")->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    int num_cells_total = 0;
    S->GetMesh("domain")->getComm()->SumAll(&num_cells, &num_cells_total, 1);
    CHECK_EQUAL(27, num_cells_total);

    int faces_middle = S->GetMesh("domain")->getSetSize(
      "middle", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
    CHECK_EQUAL(0, faces_middle);

    int cells_all =
      S->GetMesh("domain")->getSetSize("all", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    int cells_all_total = 0;
    S->GetMesh("domain")->getComm()->SumAll(&cells_all, &cells_all_total, 1);
    CHECK_EQUAL(27, cells_all_total);
  }


  TEST_FIXTURE(obs_test, ObservePoint)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "constant");
    obs_list.set<std::string>("region", "middle");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "point");

    Teuchos::ParameterList vlist;
    vlist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    auto vo = Teuchos::rcp(new VerboseObject("Test", vlist));

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    WriteStateStatistics(*S, *vo);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(2.0, observation[0], 1.e-10);
  }


  TEST_FIXTURE(obs_test, ObserveIntensiveIntegral)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "constant");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "integral");

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(2.0 * (2 * 2 * 2), observation[0], 1.e-10); // total volume
  }


  TEST_FIXTURE(obs_test, ObserveExtensiveIntegral)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "constant");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "extensive integral");

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(2.0 * 27, observation[0], 1.e-10); // number of cells
  }


  TEST_FIXTURE(obs_test, ObserveAverage)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "constant");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "average");

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(2.0, observation[0], 1.e-10); // vol weighted average
  }

  TEST_FIXTURE(obs_test, ObserveMin)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "id");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "minimum");

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(0.0, observation[0], 1.e-10); // smallest id is 0
  }

  TEST_FIXTURE(obs_test, ObserveMax)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "id");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "maximum");

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);

    CHECK_CLOSE(26, observation[0], 1.e-10); // biggest id is num cells - 1
  }


  TEST_FIXTURE(obs_test, Face)
  {
    setup();
    // integrate an observable
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "flux");
    obs_list.set<std::string>("region", "top face");
    obs_list.set<std::string>("location name", "face");
    obs_list.set<std::string>("functional", "extensive integral");
    obs_list.set("direction normalized flux", true);

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(4.0, observation[0], 1.e-10); // flux is given by face area
  }


  TEST_FIXTURE(obs_test, Face_NORMALIZED_REL_VOLUME)
  {
    setup();
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
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(1, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(4.0, observation[0], 1.e-10); // flux is given by face area
  }


  TEST_FIXTURE(obs_test, MULTI_DOF_OBS_ALL)
  {
    setup();
    // direction nomralized flux relative to region allows normalizing in an
    // outward-normal relative to a volumetric region.
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "multi_dof");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "average");

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(3, obs.get_num_vectors());

    std::vector<double> observation(3, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(0.0, observation[0], 1.e-10);
    CHECK_CLOSE(1.0, observation[1], 1.e-10);
    CHECK_CLOSE(2.0, observation[2], 1.e-10);
  }


  TEST_FIXTURE(obs_test, MULTI_DOF_OBS_ONE)
  {
    setup();
    // direction nomralized flux relative to region allows normalizing in an
    // outward-normal relative to a volumetric region.
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("variable", "multi_dof");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("functional", "average");
    obs_list.set<int>("degree of freedom", 2);

    Observable obs(obs_list);
    obs.set_comm(S->GetMesh("domain")->getComm());
    obs.Setup(S.ptr());
    obs.FinalizeStructure(S.ptr());
    CHECK_EQUAL(1, obs.get_num_vectors());

    std::vector<double> observation(3, Observable::nan);
    obs.Update(S.ptr(), observation, 0);
    CHECK_CLOSE(2.0, observation[0], 1.e-10);
  }


  TEST_FIXTURE(obs_test, FileOne)
  {
    setup();
    //  one observation in a file
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("observation output filename", "obs1.dat");
    obs_list.set<Teuchos::Array<int>>("cycles", std::vector<int>{ 0, 1 });
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


  TEST_FIXTURE(obs_test, FileTwo)
  {
    setup();
    //  one observation in a file
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("observation output filename", "obs2.dat");
    obs_list.set<Teuchos::Array<int>>("cycles", std::vector<int>{ 0, 1 });

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


  TEST_FIXTURE(obs_test, TimeIntegrated)
  {
    setup();
    //  one observation in a file
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("observation output filename", "obs3.dat");
    obs_list.set<Teuchos::Array<double>>("times", std::vector<double>{ 0, 1 });

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


  TEST_FIXTURE(obs_test, WritesNaN)
  {
    setup();
    // integrate an observable
    //  one observation in a file
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("observation output filename", "obs4.dat");
    obs_list.set<Teuchos::Array<double>>("times", std::vector<double>{ 0, 1 });

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

  TEST_FIXTURE(obs_test, FileWithModifier)
  {
    setup();
    //  one observation in a file
    Teuchos::ParameterList obs_list("my obs");
    obs_list.set<std::string>("observation output filename", "obs5.dat");
    obs_list.set<Teuchos::Array<int>>("cycles", std::vector<int>{ 0, 1 });
    obs_list.set<std::string>("variable", "linear");
    obs_list.set<std::string>("region", "all");
    obs_list.set<std::string>("location name", "cell");
    obs_list.set<std::string>("reduction", "average");

    Teuchos::ParameterList& func_list = obs_list.sublist("modifier").sublist("function-linear");
    Teuchos::Array<double> grad(1);
    grad[0] = 2.0;
    func_list.set("gradient", grad);
    func_list.set("y0", (double)0);

    {
      UnstructuredObservations obs(obs_list);
      obs.Setup(S.ptr());
      obs.MakeObservations(S.ptr());
      advance(1.0);
      obs.MakeObservations(S.ptr());
    }

    // times: 0, 1
    // values: 0, .1
    CHECK(compareFiles("obs5.dat", "test/obs5.dat.gold"));
  }

  TEST_FIXTURE(obs_domain_set_test, ObsDomainSet)
  {
    setup_domain_set();

    // must be able to deal with column observations off process, and still write
    // only on the local process if all are local
    Teuchos::ParameterList obs_list1("my obs ds1");
    obs_list1.set<std::string>("observation output filename", "obs_ds1.dat");
    obs_list1.set<Teuchos::Array<double>>("times", std::vector<double>{ 0, 1 });

    auto& obsA_list = obs_list1.sublist("observed quantities").sublist("obsA");
    obsA_list.set<std::string>("variable", "column:0-variable");
    obsA_list.set<std::string>("region", "all");
    obsA_list.set<std::string>("location name", "cell");
    obsA_list.set<std::string>("functional", "average");

    auto& obsB_list = obs_list1.sublist("observed quantities").sublist("obsB");
    obsB_list.set<std::string>("variable", "column:8-variable");
    obsB_list.set<std::string>("region", "all");
    obsB_list.set<std::string>("location name", "cell");
    obsB_list.set<std::string>("functional", "average");

    {
      UnstructuredObservations obs(obs_list1);
      obs.Setup(S.ptr());
      obs.MakeObservations(S.ptr());
      advance(0.5);
      obs.MakeObservations(S.ptr());
      advance(0.5);
      obs.MakeObservations(S.ptr());
    }

    // times: 0, 1
    // valuesA: 0,0
    // valuesB: 8,8
    CHECK(compareFiles("obs_ds1.dat", "test/obs_ds1.dat.gold"));
  }


  TEST_FIXTURE(obs_domain_set_test, ObsDomainSetSubCommunicator)
  {
    setup_domain_set();

    // observation that exists and is written on rank 1
    Teuchos::ParameterList obs_list("my obs ds2");
    obs_list.set<std::string>("observation output filename", "obs_ds2.dat");
    obs_list.set<Teuchos::Array<double>>("times", std::vector<double>{ 0, 1 });
    obs_list.set<std::string>("domain", "column:8");

    auto& obsB_list = obs_list.sublist("observed quantities").sublist("obsB");
    obsB_list.set<std::string>("variable", "column:8-variable");
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
    // valuesB: 8,8
    CHECK(compareFiles("obs_ds2.dat", "test/obs_ds2.dat.gold"));
  }
}
