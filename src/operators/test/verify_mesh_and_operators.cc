/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

// TPLs
#include "UnitTest++.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"

// Amanzi
#include "bilinear_form_registration.hh"

#include "DiffusionFixture.hh"
#include "DiffusionFixtureTests.hh"

static std::vector<std::string> argv_copy;

int main(int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << "Usage: verify_mesh_and_operators  pc_type  mesh_type  mesh_size|mesh_file  scheme  tol\n"
                 "  (req) pc_type = identity|diagonal|ifpack2: ILUT|Hypre: AMG\n"
                 "  (req) mesh_type = structured2d|structured3d|unstructured2d|unstructured3d\n"
                 "  (opt) mesh_size = positive integer\n"
                 "  (opt) mesh_file = file containing mesh\n"
                 "  (opt) scheme    = mfd|fv\n"
                 "  (opt) tol       = 1e-10\n";
    return 1;
  }
  for (int i = 1; i < argc; ++i) argv_copy.push_back(argv[i]);

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return UnitTest::RunAllTests();
}


TEST(Verify_Mesh_and_Operators) {

  std::string pc_type = argv_copy[0];
  std::string mesh_type = argv_copy[1];
  std::string mesh_file("structured2d");

  // define mesh
  int d(2), n(10);

  int argc = argv_copy.size();
  if (argc > 2) {
    if (mesh_type == "structured2d") {
      d = 2;
      n = std::stoi(argv_copy[2]);
    } else if (mesh_type == "structured3d") {
      d = 3;
      n = std::stoi(argv_copy[2]);
    } else if (mesh_type == "unstructured2d") {
      d = 2;
      mesh_file = argv_copy[2];
    } else if (mesh_type == "unstructured3d") {
      d = 3;
      mesh_file = argv_copy[2];
    } else {
      std::cout << "Incorrect #2, available options: integer | exodus file name\n";
      exit(1);
    }
  }

  // define numerical scheme
  std::string scheme("mixed");
  if (argc > 3) {
    scheme = argv_copy[3];
    if (scheme == "mfd") scheme = "mixed";
  }

  double tol(1e-10);
  if (argc > 4) {
    tol = std::stod(argv_copy[4]);
  }

  // little_k
  AmanziMesh::Entity_kind scalar_coef(AmanziMesh::Entity_kind::UNKNOWN);

  // other parameters
  bool symmetric(false);
  int order(1);
  std::string ana("00");

  // create parameter list
  Teuchos::Array<std::string> dofs(2);
  dofs[0] = "face";
  dofs[1] = "cell";

  auto plist = Teuchos::rcp(new Teuchos::ParameterList());
  plist->sublist("PK operator").sublist("mixed")
      .set<std::string>("discretization primary", "mfd: default")
      .set<Teuchos::Array<std::string> >("schema", dofs)
      .set<Teuchos::Array<std::string> >("preconditioner schema", dofs)
      .set<std::string>("nonlinear coefficient", "none");

  plist->sublist("solvers").sublist("AztecOO CG")
      .set<std::string>("iterative method", "pcg").sublist("pcg parameters")
      .set<int>("maximum number of iterations", 1000)
      .set<double>("error tolerance", 1e-12);

  plist->sublist("solvers").sublist("GMRES")
      .set<std::string>("iterative method", "gmres").sublist("gmres parameters")
      .set<int>("maximum number of iterations", 1000)
      .set<double>("error tolerance", 1e-12);

  plist->sublist("preconditioners").sublist("identity")
      .set<std::string>("preconditioning method", "identity");

  plist->sublist("preconditioners").sublist("diagonal")
      .set<std::string>("preconditioning method", "diagonal").sublist("diagonal parameters");

  plist->sublist("preconditioners").sublist("ifpack2: ILUT")
      .set<std::string>("preconditioning method", "ifpack2: ILUT").sublist("ifpack2: ILUT parameters")
      .set<double>("fact: ilut level-of-fill", 10.0)
      .set<double>("fact: drop tolerance", 0.0);

  plist->sublist("preconditioners").sublist("Hypre: AMG")
      .set<std::string>("preconditioning method", "boomer amg").sublist("boomer amg parameters")
      .set<int>("cycle applications", 2)
      .set<int>("smoother sweeps", 3)
      .set<double>("strong threshold", 0.5)
      .set<double>("tolerance", 0.0)
      .set<int>("relaxation type", 6);

  test(pc_type, "Dirichlet", mesh_file, d, n,
       scheme, symmetric, scalar_coef, tol, order, ana, 1, plist);
}

