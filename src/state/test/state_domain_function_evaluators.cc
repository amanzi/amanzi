/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorIndependentPatchFunction.hh"
#include "EvaluatorModelPatch.hh"
#include "EvaluatorAggregateBCs.hh"
#include "BCs.hh"
#include "dag_models.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/*
  Tests independent and secondary variable evaluators based on functions.
*/

SUITE(DOMAIN_FUNCTIONS)
{
  //
  // Test how sources might use MeshFunctions on patches
  //
  TEST(DOMAIN_FUNCTIONS_SOURCES)
  {
    auto comm = Amanzi::getDefaultComm();

    // Create the geometric model
    Teuchos::ParameterList regions;
    std::vector<double> low = { 0., 0., 0. };
    std::vector<double> high = { 4., 4., 4. };
    regions.sublist("left")
      .sublist("region: box")
      .set<Teuchos::Array<double>>("low coordinate", low)
      .set<Teuchos::Array<double>>("high coordinate", std::vector<double>{ 2., 4., 4. });
    regions.sublist("right")
      .sublist("region: box")
      .set<Teuchos::Array<double>>("low coordinate", std::vector<double>{ 2., 0., 0. })
      .set<Teuchos::Array<double>>("high coordinate", high);
    regions.sublist("point")
      .sublist("region: point")
      .set<Teuchos::Array<double>>("coordinate", std::vector<double>{ 1., 1., 1. });
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));

    MeshFactory meshfac(comm, gm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(mesh);

    {
      // this patch is a function-provided patch
      auto& ps = S->Require<MultiPatch<double>, MultiPatchSpace>(
        "source_left", Tags::DEFAULT, "source_left");
      ps.set_mesh(mesh);
      ps.set_entity_kind(AmanziMesh::Entity_kind::CELL);

      // NOTE: this one, we do not need to addPatch -- the region is set by the
      // function spec/MeshFunction which gets read out of the parameter list.
      //  ps.addPatch("left", AmanziMesh::Entity_kind::CELL, 1);
      auto plist = Teuchos::rcp(new Teuchos::ParameterList("source_left"));
      plist->set<std::string>("function inner list name", "well pressure [m]");

      Teuchos::ParameterList& f1 = plist->sublist("function").sublist("source left");
      f1.set<std::string>("region", "left");
      f1.sublist("well pressure [m]").set<std::string>("function type", "constant");
      f1.sublist("well pressure [m]").set<double>("value", 3.14);

      auto ps_eval = Teuchos::rcp(new EvaluatorIndependentPatchFunction(plist));
      S->SetEvaluator("source_left", Tags::DEFAULT, ps_eval);
    }

    {
      // require the dependency, called "G" in the DModel used here.
      // here G is a vector, e.g. head, on the mesh.
      S->Require<CompositeVector, CompositeVectorSpace>("G", Tags::DEFAULT, "G")
        .SetMesh(mesh)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

      auto glist = Teuchos::rcp(new Teuchos::ParameterList("G"));
      Teuchos::ParameterList& g1 = glist->sublist("function").sublist("right").sublist("function");
      g1.set<std::string>("function type", "constant");
      g1.set<double>("value", 2.0);
      glist->sublist("verbose object").set<std::string>("verbosity level", "extreme");
      auto Geval = Teuchos::rcp(new EvaluatorIndependentFunction(glist));
      S->SetEvaluator("G", Tags::DEFAULT, Geval);

      // this patch is a user-provided evaluator.  Think, for instance, a sink
      // that is a function of head.  Here we will use DModel to provide that function (2 * G)
      auto& ps = S->Require<MultiPatch<double>, MultiPatchSpace>("D", Tags::DEFAULT, "D");
      ps.set_mesh(mesh);
      //  ps.set_entity_kind(AmanziMesh::Entity_kind::CELL);

      // NOTE: here, we _DO_ have to set the region, because no one else reads a
      // parameterlist...
      ps.addPatch("right", AmanziMesh::Entity_kind::CELL, 1);

      auto dlist = Teuchos::rcp(new Teuchos::ParameterList("D"));
      dlist->set<std::string>("tag", "");
      dlist->sublist("verbose object").set<std::string>("verbosity level", "extreme");
      auto Deval = Teuchos::rcp(new EvaluatorModelPatch<DModelAccessor>(dlist));
      S->SetEvaluator("D", Tags::DEFAULT, Deval);
    }

    // Setup fields and marked as initialized
    S->Setup();
    S->Initialize();

    // set a time
    S->set_time(1.1);
    S->GetEvaluator("D", Tags::DEFAULT).Update(*S, "main");

    // what to do on HOST/DEVICE?  this may fail on host since data is on DEVICE?
    auto on_host = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace(),
      S->Get<MultiPatch<double>>("source_left", Tags::DEFAULT)[0].data);
    CHECK_CLOSE(3.14, on_host(0, 0), 1.e-10);

    std::cout << "Dptr = " << S->GetPtr<MultiPatch<double>>("D", Tags::DEFAULT) << std::endl;
    auto on_host2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
            S->Get<MultiPatch<double>>("D", Tags::DEFAULT)[0].data);
    CHECK_CLOSE(4, on_host2(0, 0), 1.e-10);
  }


  //
  // Test how BCs might use MeshFunctions on patches
  //
  TEST(DOMAIN_FUNCTIONS_BCS)
  {
    auto comm = Amanzi::getDefaultComm();

    // Create the geometric model
    Teuchos::ParameterList regions;
    std::vector<double> low = { 0., 0., 0. };
    std::vector<double> high = { 4., 4., 4. };
    regions.sublist("left")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("point", low)
      .set<Teuchos::Array<double>>("normal", std::vector<double>{ -1, 0, 0 });
    regions.sublist("right")
      .sublist("region: plane")
      .set<Teuchos::Array<double>>("point", high)
      .set<Teuchos::Array<double>>("normal", std::vector<double>{ 1, 0, 0 });
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));

    MeshFactory meshfac(comm, gm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(mesh);

    {
      // this patch is a function-provided patch
      auto& ps =
        S->Require<MultiPatch<double>, MultiPatchSpace>("bc_left", Tags::DEFAULT, "bc_left");
      // ps.set_mesh(mesh);
      // ps.set_entity_kind(AmanziMesh::Entity_kind::FACE); // these two set by BCs?
      ps.set_flag(Operators::OPERATOR_BC_DIRICHLET);
      // ps.addPatch("left", AmanziMesh::Entity_kind::FACE, 1); // this set by func_?

      auto plist = Teuchos::rcp(new Teuchos::ParameterList("bc_left"));
      plist->set<std::string>("function list name", "pressure");
      plist->set<std::string>("function inner list name", "pressure [Pa]");

      Teuchos::ParameterList& f1 = plist->sublist("pressure").sublist("source left");
      f1.set<std::string>("region", "left");
      f1.sublist("pressure [Pa]").set<std::string>("function type", "constant");
      f1.sublist("pressure [Pa]").set<double>("value", 3.14);

      auto ps_eval = Teuchos::rcp(new EvaluatorIndependentPatchFunction(plist));
      S->SetEvaluator("bc_left", Tags::DEFAULT, ps_eval);
    }

    {
      // Imagine this is a critical depth condition, where q is a function of h
      S->Require<CompositeVector, CompositeVectorSpace>("G", Tags::DEFAULT, "G")
        .SetMesh(mesh)
        ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

      // G here plays the role of h
      auto glist = Teuchos::rcp(new Teuchos::ParameterList("G"));
      Teuchos::ParameterList& g1 = glist->sublist("function").sublist("right").sublist("function");
      g1.set<std::string>("function type", "constant")
        .set<std::string>("component", "face")
        .set<double>("value", 2.0);
      glist->sublist("verbose object").set<std::string>("verbosity level", "extreme");
      auto Geval = Teuchos::rcp(new EvaluatorIndependentFunction(glist));
      S->SetEvaluator("G", Tags::DEFAULT, Geval);

      // D plays the role of q(h), computed via the DModel as 2*G
      auto& ps = S->Require<MultiPatch<double>, MultiPatchSpace>("D", Tags::DEFAULT, "D");
      // ps.set_mesh(mesh);
      // ps.set_entity_kind(AmanziMesh::Entity_kind::FACE); // these two set by BCs?
      ps.set_flag(Operators::OPERATOR_BC_NEUMANN); // but we must set the flag
      // ps.addPatch("right", AmanziMesh::Entity_kind::FACE, 1);
      ps.addPatch("right", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

      auto dlist = Teuchos::rcp(new Teuchos::ParameterList("D")); // requires no list?
      dlist->set<std::string>("tag", "");
      dlist->set<int>("flag", Operators::OPERATOR_BC_NEUMANN);
      dlist->sublist("verbose object").set<std::string>("verbosity level", "extreme");
      auto Deval = Teuchos::rcp(new EvaluatorModelPatch<DModelAccessor>(dlist));
      S->SetEvaluator("D", Tags::DEFAULT, Deval);
    }

    {
      // lastly require the BCs object that aggregates these for use by PDE_Diffusion
      auto& bc_fac = S->Require<Operators::BCs, Operators::BCs_Factory>(
        "diffusion_bcs", Tags::DEFAULT, "diffusion_bcs");
      bc_fac.set_mesh(mesh);
      bc_fac.set_entity_kind(AmanziMesh::Entity_kind::FACE);

      // and the evaluator
      auto bclist = Teuchos::rcp(new Teuchos::ParameterList("diffusion_bcs"));
      bclist->set<std::string>("tag", "");
      bclist->set<Teuchos::Array<std::string>>("dependencies",
                                               std::vector<std::string>{ "bc_left", "D" });
      bclist->sublist("verbose object").set<std::string>("verbosity level", "extreme");
      auto bceval = Teuchos::rcp(new EvaluatorAggregateBCs(bclist));
      S->SetEvaluator("diffusion_bcs", Tags::DEFAULT, bceval);
    }

    // Setup fields and marked as initialized
    S->Setup();
    S->Initialize();

    // set a time
    S->set_time(1.1);
    S->GetEvaluator("diffusion_bcs", Tags::DEFAULT).Update(*S, "main");

    // compare to make sure we get flags and values right...
    {
      const auto& bc_markers =
        S->Get<Operators::BCs>("diffusion_bcs", Tags::DEFAULT).bc_model<MemSpace_kind::HOST>();
      const auto& bc_values =
        S->Get<Operators::BCs>("diffusion_bcs", Tags::DEFAULT).bc_value<MemSpace_kind::HOST>();

      int nfaces =
        mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
      CHECK_EQUAL(nfaces, bc_markers.extent(0));
      CHECK_EQUAL(nfaces, bc_values.extent(0));

      for (int f = 0; f != nfaces; ++f) {
        auto fc = mesh->getFaceCentroid(f);
        if (std::abs(fc[0]) < 1.e-10) {
          // should be a Dirichlet face, value of 3.14
          CHECK_EQUAL(Operators::OPERATOR_BC_DIRICHLET, bc_markers(f));
          CHECK_CLOSE(3.14, bc_values(f), 1.e-10);
        } else if (std::abs(fc[0] - 4.0) < 1.e-10) {
          // shoudl be Neumann, value of 4
          CHECK_EQUAL(Operators::OPERATOR_BC_NEUMANN, bc_markers(f));
          CHECK_CLOSE(4.0, bc_values(f), 1.e-10);
        } else if (std::abs(fc[1]) < 1.e-10 || std::abs(fc[1] - 4.0) < 1.e-10 ||
                   std::abs(fc[2]) < 1.e-10 || std::abs(fc[2] - 4.0) < 1.e-10) {
          // boundary face on y or z, Neumann no flux
          CHECK_EQUAL(Operators::OPERATOR_BC_NEUMANN, bc_markers(f));
          CHECK_CLOSE(0.0, bc_values(f), 1.e-10);
        } else {
          // NO BC
          CHECK_EQUAL(Operators::OPERATOR_BC_NONE, bc_markers(f));
          CHECK_CLOSE(0.0, bc_values(f), 1.e-10);
        }
      }
    }
  }
}
