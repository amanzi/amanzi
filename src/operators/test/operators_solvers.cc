/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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
#include "AmanziTypes.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
#include "Analytic02.hh"
#include "Analytic03b.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFVwithGravity.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_DiffusionMFDwithGravity.hh"

#include "DiffusionFixtureOld.hh"

#include "cuda_decl.h"


struct TestHarness {};

template<class PDE_Diffusion_type>
void test(const Teuchos::RCP<AnalyticBase>& ana,
          const std::string& pc_type,
          const std::string& solver_type,
          const std::string& bc_type,
          const std::string& mesh_type,
          const std::string& disc_type,
          AmanziMesh::Entity_kind scalar_coef,
          double tol,
          int niters=1) {
  nvtxRangePushA("Generating mesh");
  auto plist = Teuchos::getParametersFromXmlFile("test/operator_solvers.xml");
  auto comm = Amanzi::getDefaultComm();
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(ana->dimension(),
            plist->sublist("regions"), *comm));
  AmanziMesh::MeshFactory meshfactory(comm, gm);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 10, 100);
  nvtxRangePop();

  nvtxRangePushA("Create fix");
  DiffusionFixture fix(ana, mesh_type, mesh);
  nvtxRangePop();

  std::cout << std::endl << std::endl << std::endl
            << "================================================================================" << std::endl
            << "Diffusion Test (np=" << fix.comm->getSize() << "): " << std::endl
            << "  problem name: " << ana->name() << std::endl
            << "  discretization: " << disc_type << std::endl
            << "  mesh: " << mesh_type << std::endl
            << "  preconditioner: " << pc_type << std::endl
            << "  solver: " << solver_type << std::endl
            << "--------------------------------------------------------------------------------"
            << std::endl;
  nvtxRangePushA("Discretize fix");
  fix.discretize<PDE_Diffusion_type, AmanziMesh::FACE>(disc_type);
  nvtxRangePop();
  if (scalar_coef != AmanziMesh::Entity_kind::UNKNOWN) fix.scalarCoefficient(scalar_coef);
  if (bc_type == "Dirichlet") {
    fix.setBCsDirichlet();
  } else if (bc_type == "DirichletNeumannBox") {
    fix.setBCsDirichletNeumannBox();
  } else {
    AMANZI_ASSERT(false);
  }
  nvtxRangePushA("Fix setup");
  fix.setup(pc_type, solver_type);
  nvtxRangePop();
  nvtxRangePushA("Fix go");
  for (int i=0; i!=niters-1; ++i) fix.go(0.0);
  fix.go(tol);
  nvtxRangePop();
  std::cout << "=============================================================================" << std::endl;
}

SUITE(SOLVERS_EASY) {
  //
  // This suite tests the range of linear solvers we have available on a
  // simple, standard diffusion problem that looks a lot like our real
  // problems.  It uses the FV discretization (probably should make this MFD)
  // and a simple generated mesh.
  //
  // It is designed for comparing linear solvers on a realish problem.
  //
  // This is the "easy" problem which uses constant coefficients = 1 and the
  // finite volume discretization.
  // TEST() {
  //   nvtxRangePushA("");
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "", "Dirichlet", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::FACE, 1.e-12);
  //   nvtxRangePop();
  // }

  TEST(IDENTITY) {
    nvtxRangePushA("identity easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "identity", "PCG", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-12);
    nvtxRangePop();
  }

  TEST(DIAGONAL) {
    nvtxRangePushA("diagonal easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "diagonal", "PCG", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-12);
    nvtxRangePop();
  }

  TEST(IFPACK2_ILUT_CPU_EASY) {
    nvtxRangePushA("ifpack2: ilut easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "ifpack2: ILUT", "PCG", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-12);
    nvtxRangePop();
  }

  TEST(IFPACK2_SCHWARZ_ILUT_EASY) {
    nvtxRangePushA("ifpack2: additive schwarz ilut easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "ifpack2: SCHWARZ", "PCG", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-12);
    nvtxRangePop();
  }

  TEST(IFPACK2_RILUK_EASY) {
    nvtxRangePushA("ifpack2: riluk easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "ifpack2: RILUK", "PCG", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-12);
    nvtxRangePop();
  }

  // NOTE: direct solve fluxes are not as accurate?
  TEST(AMESOS2_KLU_EASY) {
    nvtxRangePushA("amesos2: klu easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "none", "Amesos2: KLU", "Dirichlet", "Generate2D",
      "fv", AmanziMesh::Entity_kind::FACE, 1.e-7);
    nvtxRangePop();
  }

  TEST(AMESOS2_BASKER_EASY) {
    nvtxRangePushA("amesos2: basker easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "none", "Amesos2: basker", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-7);
    nvtxRangePop();
  }

  TEST(AMESOS2_SUPERLU_DIST_EASY) {
    nvtxRangePushA("amesos2: superlu dist easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "none", "Amesos2: superludist", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-7);
    nvtxRangePop();
  }

  TEST(HYPRE_BOOMERAMG_EASY) {
    nvtxRangePushA("hypre: boomer amg easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "hypre: boomer amg", "GMRES", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-7);
    nvtxRangePop();
  }
  TEST(MUELU_EASY) {
    nvtxRangePushA("muelu easy");
    auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
      ana, "MueLu", "GMRES", "Dirichlet", "Generate2D",
        "fv", AmanziMesh::Entity_kind::FACE, 1.e-12);
    nvtxRangePop();
  }
}

// SUITE(SOLVERS_HARD) {

//   // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
//   TEST(Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_ILUT) {
//     auto ana = Teuchos::rcp(new Analytic03b());
//     test<Operators::PDE_DiffusionMFD>(
//         ana, "ifpack2: SCHWARZ", "Dirichlet", "Generate2D_HiRes",
//         "mixed upwind", false, AmanziMesh::Entity_kind::FACE, 2.7e-2);
//   }

// }
