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
static int MyPID;

int main(int argc, char *argv[])
{
  if (argc < 4) {
    std::cout <<
      "Usage: verify_operators  with_pc|direct  mesh_type  mesh_size|mesh_file  scheme  tol  nloops  linsolver  test_id\n\n"
      "  (req) with_pc   = identity|diagonal|ifpack: ILUT\n"
      "                    Hypre: AMG|Hypre: Euclid\n"
      "                    Trilinos: ML|Trilinos: MueLu\n"
      "  (req) direct    = Amesos1: KLU|Amesos2: Basker|Amesos2: SuperLUDist\n\n"
      "  (req) mesh_type = structured2d|structured3d|unstructured2d|unstructured3d\n"
      "  (req) mesh_size = positive integer\n"
      "  (req) mesh_file = file containing mesh\n\n"
      "  (opt) scheme    = mfd|fv|mfd_upwind  (default mfd)\n"
      "  (opt) tol       = positive double  (default 1e-10)\n"
      "  (opt) nloops    = number of iterations  (default is 1 for linear solvers)\n"
      "  (opt) libsolver = linear solver: pcg (default) or gmres\n"
      "  (opt) test_id   = id of analytic test (default 00) \n\n"
      "Examples:\n"
      "  verify_operators \"Hypre: AMG\" structured3d 10 fv 1e-10\n"
      "  verify_operators \"Amesos1: KLU\" unstructured2d mymesh.exo mfd 1e-10 1 gmres\n"
      "  verify_operators \"identity\" unstructured2d mymesh.exo mfd 1e-10 1 gmres \"03\"\n";
    return 1;
  }
  for (int i = 1; i < argc; ++i) argv_copy.push_back(argv[i]);

  argv[0] = new char[40];
  strcpy(argv[0], "--teuchos-suppress-startup-banner");
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  MyPID = mpiSession.getRank();

  return UnitTest::RunAllTests();
}


TEST(Verify_Mesh_and_Operators) {

  std::string prec_solver = argv_copy[0];
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
      mesh_file = "structured3d";
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
    if (scheme == "mfd_upwind") scheme = "mixed upwind";
  }

  double tol(1e-10);
  if (argc > 4) {
    tol = std::stod(argv_copy[4]);
  }

  double nloops(1);
  if (argc > 5) {
    nloops = std::stoi(argv_copy[5]);
  }

  std::string linsolver("pcg");
  if (argc > 6) {
    linsolver = argv_copy[6];
  }

  std::string ana("00");
  if (argc > 7) {
    ana = argv_copy[7];
  }

  // set other parameters based on input values
  // -- little_k
  AmanziMesh::Entity_kind scalar_coef(AmanziMesh::Entity_kind::UNKNOWN);
  if (scheme == "mixed upwind") scalar_coef = AmanziMesh::Entity_kind::FACE;

  // -- symmetry
  int order(1);
  bool symmetric(true);
  if (linsolver != "pcg") symmetric = false;

  // -- parameter list for discretization
  Teuchos::Array<std::string> dofs(2);
  dofs[0] = "face";
  dofs[1] = "cell";

  auto plist = Teuchos::rcp(new Teuchos::ParameterList());
  plist->sublist("PK operator").sublist("mixed")
      .set<std::string>("discretization primary", "mfd: default")
      .set<Teuchos::Array<std::string> >("schema", dofs)
      .set<Teuchos::Array<std::string> >("preconditioner schema", dofs)
      .set<std::string>("nonlinear coefficient", "none");

  plist->sublist("PK operator").sublist("mixed upwind")
      .set<std::string>("discretization primary", "mfd: optimized for sparsity")
      .set<Teuchos::Array<std::string> >("schema", dofs)
      .set<Teuchos::Array<std::string> >("preconditioner schema", dofs)
      .set<std::string>("nonlinear coefficient", "upwind: face");

  plist->sublist("PK operator").sublist("so")
      .set<std::string>("discretization primary", "mfd: support operator")
      .set<Teuchos::Array<std::string> >("schema", dofs)
      .set<Teuchos::Array<std::string> >("preconditioner schema", dofs)
      .set<std::string>("nonlinear coefficient", "none");

  plist->sublist("PK operator").sublist("fv")
      .set<std::string>("discretization primary", "fv: default")
      .set<std::string>("schema", "cell")
      .set<std::string>("preconditioner schema", "cell")
      .set<std::string>("nonlinear coefficient", "none");

  // -- parameter list for solvers
  plist->sublist("solvers").sublist("AztecOO CG")
      .set<std::string>("iterative method", "pcg").sublist("pcg parameters")
      .set<int>("maximum number of iterations", 5000)
      .set<double>("error tolerance", tol);

  plist->sublist("solvers").sublist("GMRES")
      .set<std::string>("iterative method", "gmres").sublist("gmres parameters")
      .set<int>("maximum number of iterations", 5000)
      .set<double>("error tolerance", tol)
      .set<int>("size of Krylov space", 20)
      .set<bool>("release Krylov vectors", "false");

  plist->sublist("solvers").sublist("Amesos1: KLU")
      .set<std::string>("direct method", "amesos").sublist("amesos parameters")
      .set<std::string>("solver name", "Klu")
      .set<int>("amesos version", 1);

  plist->sublist("solvers").sublist("Amesos2: Basker")
      .set<std::string>("direct method", "amesos").sublist("amesos parameters")
      .set<std::string>("solver name", "basker")
      .set<int>("amesos version", 2);

  plist->sublist("solvers").sublist("Amesos2: SuperLUDist")
      .set<std::string>("direct method", "amesos").sublist("amesos parameters")
      .set<std::string>("solver name", "Superludist")
      .set<int>("amesos version", 2);

  // preconditioners
  // -- simple
  plist->sublist("preconditioners").sublist("identity")
      .set<std::string>("preconditioning method", "identity");

  plist->sublist("preconditioners").sublist("diagonal")
      .set<std::string>("preconditioning method", "diagonal").sublist("diagonal parameters");

  // -- Hypre
  plist->sublist("preconditioners").sublist("Hypre: Euclid")
      .set<std::string>("preconditioning method", "euclid").sublist("euclid parameters")
      .set<int>("ilu(k) fill level", 10)
      .set<bool>("rescale rows", false)
      .set<double>("ilut drop tolerance", 1e-6);

  plist->sublist("preconditioners").sublist("Hypre: AMG")
      .set<std::string>("preconditioning method", "boomer amg").sublist("boomer amg parameters")
      .set<int>("cycle applications", 2)
      .set<int>("smoother sweeps", 3)
      .set<double>("strong threshold", 0.5)
      .set<double>("tolerance", 0.0)
      .set<int>("verbosity", 0)
      .set<int>("relaxation type", 6);

  // -- Trilinos
  plist->sublist("preconditioners").sublist("Trilinos: ML")
      .set<std::string>("preconditioning method", "ml").sublist("ml parameters")
      .set<int>("ML output", 0)
      .set<int>("aggregation: nodes per aggregate", 3)
      .set<double>("aggregation: damping factor", 1.33)
      .set<double>("aggregation: threshold", 0.0)
      .set<std::string>("aggregation: type", "Uncoupled")
      .set<int>("coarse: max size", 128)
      .set<std::string>("coarse: type", "Amesos-KLU")
      .set<int>("cycle applications", 2)
      .set<int>("eigen-analysis: iterations", 10)
      .set<std::string>("eigen-analysis: type", "cg")
      .set<int>("max levels", 40)
      .set<std::string>("prec type", "MGW")
      .set<double>("smoother: damping factor", 1.0)
      .set<std::string>("smoother: pre or post", "both")
      .set<int>("smoother: sweeps", 2);

  plist->sublist("preconditioners").sublist("Trilinos: MueLu")
    .set<std::string>("preconditioning method", "muelu").sublist("muelu parameters")
    .set<int>("max levels", 25)
    .set<int>("coarse: max size", 10)
    .set<std::string>("coarse: type", "SuperLU_dist")
    .set<std::string>("verbosity", "low")
    .set<std::string>("multigrid algorithm", "sa")
    .set<std::string>("smoother: type", "RELAXATION")
    .set<std::string>("aggregation: type", "uncoupled")
    .set<int>("aggregation: min agg size", 3)
    .set<int>("aggregation: max agg size", 9).sublist("smoother: params")
      .set<std::string>("relaxation: type", "symmetric Gauss-Seidel")
      .set<int>("relaxation: sweeps", 1)
      .set<double>("relaxation: damping factor", 0.9);

  // -- ILU
  plist->sublist("preconditioners").sublist("ifpack: ILUT")
      .set<std::string>("preconditioning method", "ifpack: ILUT").sublist("ifpack: ILUT parameters")
      .set<double>("fact: ilut level-of-fill", 10.0)
      .set<double>("fact: drop tolerance", 0.0);

  plist->sublist("preconditioners").sublist("ifpack2: ILUT")
      .set<std::string>("preconditioning method", "ifpack2: ILUT").sublist("ifpack2: ILUT parameters")
      .set<double>("fact: ilut level-of-fill", 10.0)
      .set<double>("fact: drop tolerance", 0.0);

  // -- RILUK: a modified variant of the ILU(k) factorization
  // -- RILUK: a modified variant of the ILU(k) factorization
  /*
  plist->sublist("preconditioners").sublist("ifpack2: RILUK")
      .set<std::string>("preconditioning method", "ifpack2: RILUK").sublist("ifpack2: RILUK parameters")
      .set<int>("fact: iluk level-of-fill", 10)
      .set<double>("fact: drop tolerance", 0.0)
      .set<std::string>("fact: type", "KSPILUK");
  */

  // summary of options
  if (MyPID == 0) {
    std::cout << "================================================================================\n";
    if (prec_solver.substr(0, 6) == "Amesos")
      plist->sublist("solvers").sublist(prec_solver).print(std::cout, 0, true, false);
    else
      plist->sublist("preconditioners").sublist(prec_solver).print(std::cout, 0, true, false);
  }

  double error_tol(10 * tol); 
  if (ana == "03") error_tol = 0.1;

  test(prec_solver, "Dirichlet", mesh_file, d, n,
       scheme, symmetric, scalar_coef, error_tol, order, ana, nloops, plist);
}

