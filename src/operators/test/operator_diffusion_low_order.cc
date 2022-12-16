/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
#include "Analytic02.hh"
#include "Analytic03b.hh"

#include "OperatorDefs.hh"
// #include "UpwindSecondOrder.hh"
// #include "Verification.hh"

#include "DiffusionFixture.hh"
#include "DiffusionFixtureTests.hh"

#define FV 1
#define SO 1
#define MFD 1
#define NLFV 1
#define ASSEMBLING 1


/* *****************************************************************
* Exactness test for mixed diffusion solver.
*/
SUITE(DIFFUSION)
{
  // NOTE: we attempt combinations of:
  //    discretization (FV, MFD, NLFV, NLFV_BndFaces)
  //    boundary condition type (Dirichlet, DirichletNeumann)
  //    mesh:
  //            2D structured
  //            3D structured
  //            random10.exo (a square with wiggled nodes),
  //            median7x8_filtered.exo (square with polyhedra, mostly hexes, some others)
  //    preconditioner (identity, diagonal, ILUT)
  //    coefficient:
  //            tensor = 1, scalar = None
  //            tensor = 1, scalar != 1 (FACE)
  //            tensor = 1, scalar != 1 (CELL)
  //            tensor !=1, scalar = None
  //    analytic problem:
  //            Analytic00_Linear: coef=1, linear solution (exact for all methods!)
  //            Analytic00_Quadratic: coef=1, quadratic solution (exact for second order methods!)

  //            Analytic00b: 3D variant of Analytic00
  //            Analytic03: discontinuous scalar Tensor coef, no scalar coef
  //            Analytic03b: same as 03, but with Tensor coef=1 and scalar coef as the scalar
  //

  // Not all combinations make sense, obviously.  2D problems must be run with
  // 2D meshes, and 3D with 3D.
#if FV
  TEST(EXACT)
  {
    std::cout << "Test: "
              << "EXACT" << std::endl;
    DiffusionFixture fix(Teuchos::null);
    fix.Init(2, 10, "structured1d");

    std::cout << std::endl
              << "Diffusion Forward Apply Test (np=" << fix.comm->NumProc() << "): "
              << "fv"
              << ", "
              << "structured1d" << std::endl
              << "--------------------------------------------------------------------------------"
              << std::endl;

    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, 1));

    fix.Discretize("fv", AmanziMesh::FACE);
    fix.SetBCsDirichlet();
    fix.Setup("diagonal", true);

    fix.global_op->Init();
    fix.op->UpdateMatrices(Teuchos::null, fix.solution.ptr());

    CompositeVector& rhs = *fix.global_op->rhs();
    {
      int ncells =
        fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& rhs_c = *rhs.ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        rhs_c[0][c] += fix.ana->source_exact(xc, 0.0) * fix.mesh->cell_volume(c);
      }
    }

    fix.op->ApplyBCs(true, true, true);

    {
      int ncells =
        fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& soln_c = *fix.solution->ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        soln_c[0][c] = fix.ana->pressure_exact(xc, 0.0);
      }
    }

    CompositeVector res(fix.global_op->DomainMap());
    fix.global_op->ComputeResidual(*fix.solution, res);

    double norminf;
    res.NormInf(&norminf);
    std::cout << "Residual norm = " << norminf << std::endl;
    CHECK_CLOSE(0.0, norminf, 1e-10);
  }
#endif

#if MFD
  TEST(EXACT_MFD)
  {
    std::cout << "Test: "
              << "EXACT_MFD" << std::endl;
    DiffusionFixture fix(Teuchos::null);
    fix.Init(2, 10, "structured1d");

    std::cout << std::endl
              << "================================================================================"
              << std::endl
              << "Diffusion Forward Apply Test (np=" << fix.comm->NumProc() << "): "
              << "MFD"
              << ", "
              << "structured1d" << std::endl
              << "--------------------------------------------------------------------------------"
              << std::endl;

    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, 1));
    fix.Discretize("mixed upwind", AmanziMesh::FACE);
    fix.SetBCsDirichlet();
    fix.Setup("diagonal", true);

    fix.global_op->Init();
    fix.op->UpdateMatrices(Teuchos::null, fix.solution.ptr());

    CompositeVector& rhs = *fix.global_op->rhs();
    {
      int ncells =
        fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& rhs_c = *rhs.ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        rhs_c[0][c] += fix.ana->source_exact(xc, 0.0) * fix.mesh->cell_volume(c);
      }
    }

    fix.op->ApplyBCs(true, true, true);

    {
      int ncells =
        fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& soln_c = *fix.solution->ViewComponent("cell");
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        soln_c[0][c] = fix.ana->pressure_exact(xc, 0.0);
      }
    }

    {
      int nfaces =
        fix.mesh->num_entities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      auto& soln_f = *fix.solution->ViewComponent("face", false);
      for (int f = 0; f != nfaces; ++f) {
        const auto& xf = fix.mesh->face_centroid(f);
        soln_f[0][f] = fix.ana->pressure_exact(xf, 0.0);
      }
    }

    CompositeVector res(fix.global_op->DomainMap());
    fix.global_op->ComputeResidual(*fix.solution, res);
    double norminf;
    res.NormInf(&norminf);
    std::cout << "Residual norm = " << norminf << std::endl;
    CHECK_CLOSE(0.0, norminf, 1e-10);
  }
#endif


  //
  // Analytic00_Linear: tests exactness for linear
  // polynomial with coefficient=1
  //
  // Exact solution
#if FV
  TEST(Analytic00_Linear1_FV_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_FV_Dirichlet_structured2d_identity" << std::endl;
    test("identity",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1.e-12);
  }
  TEST(Analytic00_LinearGravity1_FV_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_FV_Dirichlet_structured2d_identity" << std::endl;
    testWGravity(1.1,
                 "identity",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1.e-12);
  }
  TEST(Analytic00_Linear1_FV_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_FV_Dirichlet_structured2d_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1.e-12);
  }

  TEST(Analytic00_LinearGravity1_FV_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_FV_Dirichlet_structured2d_diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1.e-12);
  }
#endif
#if SO
  TEST(Analytic00_Linear1_SO_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_SO_Dirichlet_structured2d_identity" << std::endl;
    test("identity",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "so",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1.e-12);
  }
#endif
#if MFD
  TEST(Analytic00_Linear1_MFD_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_MFD_Dirichlet_structured2d_identity" << std::endl;
    test("identity",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_MFD_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_MFD_Dirichlet_structured2d_identity" << std::endl;
    testWGravity(1.1,
                 "identity",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_Linear1_MFD_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_MFD_Dirichlet_structured2d_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_MFD_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_MFD_Dirichlet_structured2d_diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN);
  }
#endif
#if NLFV
  TEST(Analytic00_Linear1_NLFV_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFV_Dirichlet_structured2d_identity" << std::endl;
    test("identity",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_NLFV_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFV_Dirichlet_structured2d_identity" << std::endl;
    testWGravity(1.1,
                 "identity",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12);
  }
  TEST(Analytic00_Linear1_NLFVBFace_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFVBFace_Dirichlet_structured2d_identity" << std::endl;
    test("identity",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-11);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_Dirichlet_structured2d_identity)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFVBFace_Dirichlet_structured2d_identity" << std::endl;
    testWGravity(1.1,
                 "identity",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-11);
  }
  TEST(Analytic00_Linear1_NLFV_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFV_Dirichlet_structured2d_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12);
  }
  TEST(Analytic00_LinearGravity1_NLFV_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFV_Dirichlet_structured2d_diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12);
  }
  TEST(Analytic00_Linear1_NLFVBFace_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFVBFace_Dirichlet_structured2d_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-10);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFVBFace_Dirichlet_structured2d_diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-10);
  }
#endif


  //
  // change PC to an assembling PC, change BCs
  //
  // Exact solution
#if ASSEMBLING

#  if FV
  TEST(Analytic00_Linear1_FV_Dirichlet_structured2d_ILU)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_FV_Dirichlet_structured2d_ILU" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12);
  }
  TEST(Analytic00_LinearGravity1_FV_Dirichlet_structured2d_ILU)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_FV_Dirichlet_structured2d_ILU" << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12);
  }
#  endif
#  if MFD
  TEST(Analytic00_Linear1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN);
  }
#  endif
#  if NLFV
  TEST(Analytic00_Linear1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12);
  }
  TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12);
  }
  TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12);
  }
#  endif
#endif

  //
  // with uniform kr
  //
  // Exact solution
#if FV
  TEST(Analytic00_Linearkr_FV_Dirichlet_structured2d_Diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linearkr_FV_Dirichlet_structured2d_Diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::FACE,
         1e-12);
  }
  TEST(Analytic00_LinearGravitykr_FV_Dirichlet_structured2d_Diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravitykr_FV_Dirichlet_structured2d_Diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::FACE,
                 1e-12);
  }
#endif
#if MFD
  TEST(Analytic00_Linearkr_MFD_Dirichlet_structured2d_Diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linearkr_MFD_Dirichlet_structured2d_Diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "mixed upwind",
         true,
         AmanziMesh::Entity_kind::FACE);
  }
  TEST(Analytic00_LinearGravitykr_MFD_Dirichlet_structured2d_Diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravitykr_MFD_Dirichlet_structured2d_Diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "mixed upwind",
                 true,
                 AmanziMesh::Entity_kind::FACE);
  }
#endif
#if NLFV
  TEST(Analytic00_Linearkr_NLFV_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linearkr_NLFV_Dirichlet_structured2d_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv",
         true,
         AmanziMesh::Entity_kind::FACE,
         1e-12);
  }
  TEST(Analytic00_LinearGravitykr_NLFV_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravitykr_NLFV_Dirichlet_structured2d_diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 true,
                 AmanziMesh::Entity_kind::FACE,
                 1e-12);
  }
  TEST(Analytic00_Linearkr_NLFVBFace_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_Linearkr_NLFVBFace_Dirichlet_structured2d_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::FACE,
         1e-10);
  }
  TEST(Analytic00_LinearGravitykr_NLFVBFace_Dirichlet_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravitykr_NLFVBFace_Dirichlet_structured2d_diagonal"
              << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "Dirichlet",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::FACE,
                 1e-10);
  }
#endif

  //
  // with uniform but non-1 tensor K
  //
  // Exact solution
#if FV
  TEST(Analytic00_LinearK_FV_DirichletNeumann_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearK_FV_DirichletNeumann_structured2d_diagonal" << std::endl;
    test("diagonal",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1.e-12);
  }
  TEST(Analytic00_LinearGravityK_FV_DirichletNeumann_structured2d_diagonal)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravityK_FV_DirichletNeumann_structured2d_diagonal" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1.e-12);
  }
#endif
#if MFD
  TEST(Analytic00_LinearK_MFD_DirichletNeumann_structured2d_diagona)
  {
    std::cout << "Test: "
              << "Analytic00_LinearK_MFD_DirichletNeumann_structured2d_diagona" << std::endl;
    test("diagonal",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravityK_MFD_DirichletNeumann_structured2d_diagona)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravityK_MFD_DirichletNeumann_structured2d_diagona" << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN);
  }
#endif
#if NLFV
  TEST(Analytic00_LinearK_NLFV_DirichletNeumann_structured2d_diagona)
  {
    std::cout << "Test: "
              << "Analytic00_LinearK_NLFV_DirichletNeumann_structured2d_diagona" << std::endl;
    test("diagonal",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12);
  }
  TEST(Analytic00_LinearGravityK_NLFV_DirichletNeumann_structured2d_diagona)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravityK_NLFV_DirichletNeumann_structured2d_diagona"
              << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12);
  }
  TEST(Analytic00_LinearK_NLFVBFace_DirichletNeumann_structured2d_diagona)
  {
    std::cout << "Test: "
              << "Analytic00_LinearK_NLFVBFace_DirichletNeumann_structured2d_diagona" << std::endl;
    test("diagonal",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         2e-10);
  }
  TEST(Analytic00_LinearGravityK_NLFVBFace_DirichletNeumann_structured2d_diagona)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravityK_NLFVBFace_DirichletNeumann_structured2d_diagona"
              << std::endl;
    testWGravity(1.1,
                 "diagonal",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 2e-10);
  }
#endif

  //
  // with wiggled mesh
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
#if ASSEMBLING
#  if FV
  TEST(Analytic00_Linear1_FV_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_FV_DirichletNeumann_Wiggled_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/random10.exo",
         2,
         -1,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         2.e-2);
  }
  TEST(Analytic00_LinearGravity1_FV_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_FV_DirichletNeumann_Wiggled_ifpack2_ILUT" << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/random10.exo",
                 2,
                 -1,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 2.e-2);
  }
#  endif
#  if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic00_Linear1_MFD_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_MFD_DirichletNeumann_Wiggled_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/random10.exo",
         2,
         -1,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         2e-2);
  }
  TEST(Analytic00_LinearGravity1_MFD_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_MFD_DirichletNeumann_Wiggled_ifpack2_ILUT" << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/random10.exo",
                 2,
                 -1,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 2e-2);
  }
#  endif
#  if NLFV
  TEST(Analytic00_Linear1_NLFV_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFV_DirichletNeumann_Wiggled_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/random10.exo",
         2,
         -1,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-11,
         1,
         "00",
         15);
  }
  TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFV_DirichletNeumann_Wiggled_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/random10.exo",
                 2,
                 -1,
                 "nlfv",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-11,
                 1,
                 "00",
                 15);
  }
  TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFVBFace_DirichletNeumann_Wiggled_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/random10.exo",
         2,
         -1,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-11,
         1,
         "00",
         15);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumann_Wiggled_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFVBFace_DirichletNeumann_Wiggled_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/random10.exo",
                 2,
                 -1,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-11,
                 1,
                 "00",
                 15);
  }
#  endif

  //
  // with polyhedral mesh
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
#  if FV
  TEST(Analytic00_Linear1_FV_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_FV_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/median7x8_filtered.exo",
         2,
         -1,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         3e-2);
  }
  TEST(Analytic00_LinearGravity1_FV_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_FV_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/median7x8_filtered.exo",
                 2,
                 -1,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 3e-2);
  }
#  endif
#  if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic00_Linear1_MFD_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_MFD_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/median7x8_filtered.exo",
         2,
         -1,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         2.7e-2);
  }
  TEST(Analytic00_LinearGravity1_MFD_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_MFD_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/median7x8_filtered.exo",
                 2,
                 -1,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 2.8e-2);
  }
#  endif
#  if NLFV
  TEST(Analytic00_Linear1_NLFV_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFV_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/median7x8_filtered.exo",
         2,
         -1,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1.e-12,
         1,
         "00",
         15);
  }
  TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFV_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/median7x8_filtered.exo",
                 2,
                 -1,
                 "nlfv",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 1,
                 "00",
                 15);
  }
  TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Linear1_NLFVBFace_DirichletNeumann_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/median7x8_filtered.exo",
         2,
         -1,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12,
         1,
         "00",
         15);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumann_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_LinearGravity1_NLFVBFace_DirichletNeumann_Poly_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/median7x8_filtered.exo",
                 2,
                 -1,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 1,
                 "00",
                 15);
  }
#  endif

  //
  // Analytic00_Quadratic: test case for quadratic equations
  // polynomial with coefficient=1
  //
#  if FV
  TEST(Analytic00_Quadratic1_FV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Quadratic1_FV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.006,
         2);
  }
  TEST(Analytic00_QuadraticGravity1_FV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_QuadraticGravity1_FV_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "fv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 0.006,
                 2);
  }
#  endif
#  if MFD
  TEST(Analytic00_Quadratic1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Quadratic1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.006,
         2);
  }
  TEST(Analytic00_QuadraticGravity1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_QuadraticGravity1_MFD_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 0.006,
                 2);
  }
#  endif
#  if NLFV
  TEST(Analytic00_Quadratic1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Quadratic1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.006,
         2);
  }
  TEST(Analytic00_QuadraticGravity1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_QuadraticGravity1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 0.006,
                 2);
  }
  TEST(Analytic00_Quadratic1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_Quadratic1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.006,
         2);
  }
  TEST(Analytic00_QuadraticGravity1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic00_QuadraticGravity1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(1.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 0.006,
                 2);
  }
#  endif

  //
  // Analytic02: Tensor diffusion.  Note FV cannot be exact for this problem.
  // This test replaces old operator_diffusion.cc mixed tests.
#  if FV
  TEST(Analytic02_FV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_FV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "fv",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.004,
         1,
         "02");
  }
#  endif

#  if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic02_MFD_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_MFD_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumannRobin",
         "structured2d",
         2,
         10,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12,
         1,
         "02");
  }
  TEST(Analytic02_Gravity_MFD_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_MFD_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 2,
                 "02");
  }

  // the original test tested on circle-quad
  TEST(Analytic02_MFD_DirichletNeumann_CircleQuad_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_MFD_DirichletNeumann_CircleQuad_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/circle_quad10.exo",
         2,
         -1,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12,
         1,
         "02");
  }
  TEST(Analytic02_Gravity_MFD_DirichletNeumann_CircleQuad_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_MFD_DirichletNeumann_CircleQuad_ifpack2_ILUT" << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/circle_quad10.exo",
                 2,
                 -1,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 1,
                 "02");
  }

  // test in 3D
  TEST(Analytic02_MFD_DirichletNeumann_structured3d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_MFD_DirichletNeumann_structured3d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured3d",
         3,
         5,
         "mixed",
         true,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12,
         2,
         "02");
  }
  TEST(Analytic02_Gravity_MFD_DirichletNeumann_structured3d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_MFD_DirichletNeumann_structured3d_ifpack2_ILUT" << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured3d",
                 3,
                 5,
                 "mixed",
                 true,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 1,
                 "02");
  }
#  endif

  //
  // With NLFV
  //
#  if NLFV
  TEST(Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12,
         1,
         "02",
         20);
  }
  TEST(Analytic02_Gravity_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 1,
                 "02",
                 20);
  }

  // the original test tested on circle-quad
  TEST(Analytic02_NLFV_DirichletNeumann_CircleQuad_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_NLFV_DirichletNeumann_CircleQuad_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/circle_quad10.exo",
         2,
         -1,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-11,
         1,
         "02",
         50);
  }
  TEST(Analytic02_Gravity_NLFV_DirichletNeumann_CircleQuad_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_NLFV_DirichletNeumann_CircleQuad_ifpack2_ILUT" << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/circle_quad10.exo",
                 2,
                 -1,
                 "nlfv",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-11,
                 1,
                 "02",
                 50);
  }

  // test in 3D
  //
  // These probably should pass?  Not sure why they don't.
  TEST(Analytic02_NLFV_DirichletNeumann_structured3d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_NLFV_DirichletNeumann_structured3d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured3d",
         3,
         5,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-11,
         1,
         "02",
         20);
  }
  TEST(Analytic02_Gravity_NLFV_DirichletNeumann_structured3d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_NLFV_DirichletNeumann_structured3d_ifpack2_ILUT" << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured3d",
                 3,
                 5,
                 "nlfv",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1.e-11,
                 1,
                 "02",
                 20);
  }

  //
  // With NLFV and Boundary Faces
  //
  TEST(Analytic02_NLFVwithBndFaces_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_NLFVwithBndFaces_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-12,
         1,
         "02",
         15);
  }
  TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumann_structured2d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumann_structured2d_ifpack2_ILUT"
              << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured2d",
                 2,
                 10,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-12,
                 1,
                 "02",
                 15);
  }

  // the original test tested on circle-quad
  TEST(Analytic02_NLFVwithBndFaces_DirichletNeumann_CircleQuad_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_NLFVwithBndFaces_DirichletNeumann_CircleQuad_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "test/circle_quad10.exo",
         2,
         -1,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1e-10,
         1,
         "02",
         50);
  }
  TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumann_CircleQuad_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumann_CircleQuad_ifpack2_ILUT"
              << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "test/circle_quad10.exo",
                 2,
                 -1,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-10,
                 1,
                 "02",
                 50);
  }

  // test in 3D
  TEST(Analytic02_NLFVwithBndFaces_DirichletNeumann_structured3d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_NLFVwithBndFaces_DirichletNeumann_structured3d_ifpack2_ILUT"
              << std::endl;
    test("ifpack2: ILUT",
         "DirichletNeumann",
         "structured3d",
         3,
         5,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         1.e-11,
         1,
         "02",
         20);
  }
  TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumann_structured3d_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumann_structured3d_ifpack2_ILUT"
              << std::endl;
    testWGravity(0.1,
                 "ifpack2: ILUT",
                 "DirichletNeumann",
                 "structured3d",
                 3,
                 5,
                 "nlfv with bfaces",
                 false,
                 AmanziMesh::Entity_kind::UNKNOWN,
                 1e-11,
                 1,
                 "02",
                 20);
  }
#  endif
#endif

  //
  // Analytic03b tests non-constant, scalar coefficients
  //
  // Note this has no gravity.
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
#if FV
  TEST(Analytic03b_Linear1_FV_Dirichlet_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_FV_Dirichlet_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         50,
         "fv",
         false,
         AmanziMesh::Entity_kind::FACE,
         2.7e-2,
         1,
         "03");
  }
#endif
#if MFD
  TEST(Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         50,
         "mixed upwind",
         false,
         AmanziMesh::Entity_kind::FACE,
         2.7e-2,
         2,
         "03");
  }
#endif
#if NLFV
  TEST(Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::FACE,
         3e-2,
         2,
         "03",
         10);
  }
  TEST(Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_ILUT)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::FACE,
         3e-2,
         2,
         "03",
         10);
  }
#endif

#if FV
  TEST(Analytic03b_Linear1_FV_Dirichlet_Poly_ifpack2_diagonal)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_FV_Dirichlet_Poly_ifpack2_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "fv",
         false,
         AmanziMesh::Entity_kind::FACE,
         2.7e-2,
         1,
         "03");
  }
#endif
#if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_diagonal)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "mixed upwind",
         false,
         AmanziMesh::Entity_kind::FACE,
         2.7e-2,
         1,
         "03");
  }
#endif
#if NLFV
  TEST(Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_diagonal)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::FACE,
         3e-2,
         2,
         "03",
         10);
  }
  TEST(Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_diagonal)
  {
    std::cout << "Test: "
              << "Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_diagonal" << std::endl;
    test("diagonal",
         "Dirichlet",
         "structured2d",
         2,
         10,
         "nlfv with bfaces",
         false,
         AmanziMesh::Entity_kind::FACE,
         3e-2,
         2,
         "03",
         10);
  }
#endif
}
