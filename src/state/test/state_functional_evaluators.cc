/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondaryMonotypeFromFunction.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/*
  Tests independent and secondary variable evaluators based on functions.
*/

SUITE(EVALS)
{
  TEST(EVALS_INDPENDENT_FROMFUNCTION)
  {
    auto comm = Amanzi::getDefaultComm();

    // Create the geometric model
    Teuchos::ParameterList regions;
    std::vector<double> low = { 0., 0., 0. };
    std::vector<double> high = { 4., 4., 4. };
    regions.sublist("ALL")
      .sublist("region: box")
      .set<Teuchos::Array<double>>("low coordinate", low)
      .set<Teuchos::Array<double>>("high coordinate", high);
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));

    MeshFactory meshfac(comm, gm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(mesh);

    // Independent variable evaluator
    // -- Field A and its evaluator
    // A = c1 + (x - x0) dot grad
    Teuchos::ParameterList A_list("fa");
    A_list.set("evaluator name", "fa");
    A_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    auto& A_func_list = A_list.sublist("function");
    auto& A_func_rlist = A_func_list.sublist("myregion");
    A_func_rlist.set("region", "ALL");
    A_func_rlist.set("component", "cell");
    auto& A_func_flist = A_func_rlist.sublist("function");
    auto& A_func_lin_list = A_func_flist.sublist("function-linear");

    std::vector<double> x0 = { 1.0, 2.0, 3.0, 4.0 };
    std::vector<double> grad = { 0.5, 1.0, 2.0, 2.0 };
    A_func_lin_list.set<double>("y0", 1.0);
    A_func_lin_list.set<Teuchos::Array<double>>("x0", x0);
    A_func_lin_list.set<Teuchos::Array<double>>("gradient", grad);

    S->Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fa_eval = Teuchos::rcp(new EvaluatorIndependentFunction(A_list));
    S->SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    // Setup fields and marked as initialized
    S->Setup();
    S->Initialize();

    // set a time
    S->set_time(1.1);

    // calculate field A
    fa_eval->Update(*S, "main");
    const Epetra_MultiVector& fa =
      *S->GetW<CompositeVector>("fa", Tags::DEFAULT, "fa").ViewComponent("cell");

    // eval point of 0th cell is:
    // {  1.1, 1, 1, 1 }
    // making the value
    // (  { 1.1, 1, 1, 1 } - {1.0, 2.0, 3.0, 4.0} ) *  {0.5, 1.0, 2.0, 2.0}  + 1
    // = (.05 -1 -4  -6 + 1) = -9.95
    CHECK_CLOSE(-9.95, fa[0][0], 1e-12);
  }


  TEST(EVALS_SECONDARY_FROMFUNCTION)
  {
    std::cout << "\n\nNEXT test\n";
    auto comm = Amanzi::getDefaultComm();

    // Create the geometric model
    Teuchos::ParameterList regions;
    std::vector<double> low = { 0., 0., 0. };
    std::vector<double> high = { 4., 4., 4. };
    regions.sublist("ALL")
      .sublist("region: box")
      .set<Teuchos::Array<double>>("low coordinate", low)
      .set<Teuchos::Array<double>>("high coordinate", high);
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));

    MeshFactory meshfac(comm, gm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(mesh);

    // Independent variable evaluator as before
    // -- Field A and its evaluator
    // A = c1 + (x - x0) dot grad
    Teuchos::ParameterList A_list("fa");
    A_list.set("evaluator name", "fa");
    A_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    auto& A_func_list = A_list.sublist("function");
    auto& A_func_rlist = A_func_list.sublist("myregion");
    A_func_rlist.set("region", "ALL");
    A_func_rlist.set("component", "cell");
    auto& A_func_flist = A_func_rlist.sublist("function");
    auto& A_func_lin_list = A_func_flist.sublist("function-linear");

    std::vector<double> x0 = { 1.0, 2.0, 3.0, 4.0 };
    std::vector<double> grad = { 0.5, 1.0, 2.0, 2.0 };
    A_func_lin_list.set<double>("y0", 1.0);
    A_func_lin_list.set<Teuchos::Array<double>>("x0", x0);
    A_func_lin_list.set<Teuchos::Array<double>>("gradient", grad);

    S->Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fa_eval = Teuchos::rcp(new EvaluatorIndependentFunction(A_list));
    S->SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    // -- Field B and its evaluator as a function of A
    Teuchos::ParameterList B_list("fb");
    B_list.set("evaluator name", "fb");
    B_list.set("tag", "");
    B_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");

    auto deps = std::vector<std::string>{ "fa" };
    B_list.set<Teuchos::Array<std::string>>("dependencies", deps)
      .set<bool>("dependency tags are my tag", true);
    auto& B_func_list = B_list.sublist("function");
    auto& B_func_llist = B_func_list.sublist("function-linear");

    auto x02 = std::vector<double>{ -1 };
    auto grad2 = std::vector<double>{ 2. };
    double y02 = 1.;
    B_func_llist.set<double>("y0", y02);
    B_func_llist.set<Teuchos::Array<double>>("x0", x02);
    B_func_llist.set<Teuchos::Array<double>>("gradient", grad2);

    S->Require<CompositeVector, CompositeVectorSpace>("fb", Tags::DEFAULT, "fb")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fb_eval = Teuchos::rcp(new EvaluatorSecondaryMonotypeFromFunction(B_list));
    S->SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // Setup fields and marked as initialized
    S->Setup();
    S->Initialize();

    // set a time
    S->set_time(1.1);

    // calculate field B
    fb_eval->Update(*S, "main");
    const auto& fb = *S->Get<CompositeVector>("fb", Tags::DEFAULT).ViewComponent("cell");

    // eval point of 0th cell is:
    // {  1.1, 1, 1, 1 }
    // making the value for A
    // (  { 1.1, 1, 1, 1 } - {1.0, 2.0, 3.0, 4.0} ) *  {0.5, 1.0, 2.0, 2.0}  + 1
    // = (.05 -1 -4  -6 + 1) = -9.95 and then the value for B = (-9.95 - -1)
    // * 2.0 + 1.0 = -16.9
    CHECK_CLOSE(-16.9, fb[0][0], 1e-12);
  }
}
