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

void test(const std::string& pc_type,
          const std::string& bc_type,
          const std::string& mesh_type,
          const std::string& disc_type,
          bool symmetric,
          AmanziMesh::Entity_kind scalar_coef, 
          double tol = 1.0e-12,
          int order = 1,
          const std::string& ana = "00",
          int niters = 1)
{
  DiffusionFixture fix(2, mesh_type);
  std::cout << std::endl << std::endl << std::endl
            << "================================================================================" << std::endl
            << "Diffusion Test (np=" << fix.comm->NumProc() << "): " << disc_type 
            << ", PC: " << pc_type << ", mesh=" << mesh_type << std::endl
            << "--------------------------------------------------------------------------------"
            << std::endl;

  if (ana == "00") 
    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, order));
  else if (ana == "02") 
    fix.ana = Teuchos::rcp(new Analytic02(fix.mesh));
  else if (ana == "03") 
    fix.ana = Teuchos::rcp(new Analytic03(fix.mesh));

  fix.Discretize(disc_type, scalar_coef);
  std::cout << "Discretize done" << std::endl;

  if (bc_type == "Dirichlet") {
    fix.SetBCsDirichlet();
  } else if (bc_type == "DirichletNeumannBox") {
    fix.SetBCsDirichletNeumannBox();
  } else {
    AMANZI_ASSERT(false);
  }      
  fix.Setup(pc_type, symmetric);
  for (int i = 0; i < niters - 1; ++i) fix.Go(0.0);
  fix.Go(tol);
  std::cout << "=============================================================================" << std::endl;
}


void testWGravity(double gravity,
                  const std::string& pc_type,
                  const std::string& bc_type,
                  const std::string& mesh_type,
                  const std::string& disc_type,
                  bool symmetric,
                  AmanziMesh::Entity_kind scalar_coef,
                  double tol = 1.0e-12,
                  int order = 1,
                  const std::string& ana = "00",
                  int niters = 1)
{
  DiffusionFixture fix(2, mesh_type);
  std::cout << std::endl << std::endl << std::endl
            << "================================================================================" << std::endl
            << "DiffusionWithGravity Test (np=" << fix.comm->NumProc() << "): "
            << disc_type << ", PC: " << pc_type << ", mesh=" << mesh_type << std::endl
            << "--------------------------------------------------------------------------------"
            << std::endl;

  if (ana == "00") 
    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, order, gravity));
  else if (ana == "02") 
    fix.ana = Teuchos::rcp(new Analytic02(fix.mesh, gravity));
  else if (ana == "03") 
    fix.ana = Teuchos::rcp(new Analytic03(fix.mesh));

  fix.DiscretizeWithGravity(disc_type, gravity, scalar_coef);

  if (bc_type == "Dirichlet") {
    fix.SetBCsDirichlet();
  } else if (bc_type == "DirichletNeumannBox") {
    fix.SetBCsDirichletNeumannBox();
  } else {
    AMANZI_ASSERT(false);
  }      
  fix.Setup(pc_type, symmetric);
  for (int i = 0; i < niters - 1; ++i) fix.Go(0.0);
  fix.Go(tol);
  std::cout << "=============================================================================" << std::endl;
}

#define FV 1
#define MFD 1
#define NLFV 1
#define ASSEMBLING 1


/* *****************************************************************
* Exactness test for mixed diffusion solver.
*/
SUITE(DIFFUSION) {
  // NOTE: we attempt combinations of:
  //    discretization (FV, MFD, NLFV, NLFV_BndFaces)
  //    boundary condition type (Dirichlet, DirichletNeumannBox)
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
  TEST(EXACT) {
    DiffusionFixture fix(2, "Generate1D");
    std::cout << std::endl
              << "Diffusion Forward Apply Test (np=" << fix.comm->NumProc() << "): " << "fv" << ", "
              << "Generate1D" << std::endl
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
      int ncells = fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& rhs_c = *rhs.ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        rhs_c[0][c] += fix.ana->source_exact(xc, 0.0) * fix.mesh->cell_volume(c);
      }
    }

    fix.op->ApplyBCs(true, true, true);

    {
      int ncells = fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
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
  TEST(EXACT_MFD) {
    DiffusionFixture fix(2, "Generate1D");
    std::cout << std::endl
              << "================================================================================" << std::endl
              << "Diffusion Forward Apply Test (np=" << fix.comm->NumProc() << "): " << "MFD" << ", "
              << "Generate1D" << std::endl
              << "--------------------------------------------------------------------------------" << std::endl;
    
    fix.ana = Teuchos::rcp(new Analytic00(fix.mesh, 1));
    fix.Discretize("mixed upwind", AmanziMesh::FACE);
    fix.SetBCsDirichlet();
    fix.Setup("diagonal", true);

    fix.global_op->Init();
    fix.op->UpdateMatrices(Teuchos::null, fix.solution.ptr());

    CompositeVector& rhs = *fix.global_op->rhs();
    {
      int ncells = fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& rhs_c = *rhs.ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        rhs_c[0][c] += fix.ana->source_exact(xc, 0.0) * fix.mesh->cell_volume(c);
      }
    }

    fix.op->ApplyBCs(true, true, true);

    {
      int ncells = fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      auto& soln_c = *fix.solution->ViewComponent("cell");
      for (int c = 0; c != ncells; ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        soln_c[0][c] = fix.ana->pressure_exact(xc, 0.0);
      }
    }

    {
      int nfaces = fix.mesh->num_entities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
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
  TEST(Analytic00_Linear1_FV_Dirichlet_Generate2D_identity) {
    test("identity", "Dirichlet", "Generate2D",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  }
  TEST(Analytic00_LinearGravity1_FV_Dirichlet_Generate2D_identity) {
    testWGravity(1.1, "identity", "Dirichlet", "Generate2D",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  }
  TEST(Analytic00_Linear1_FV_Dirichlet_Generate2D_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  }

  TEST(Analytic00_LinearGravity1_FV_Dirichlet_Generate2D_diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  } 
#endif
#if MFD  
  TEST(Analytic00_Linear1_MFD_Dirichlet_Generate2D_identity) {
    test("identity", "Dirichlet", "Generate2D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_MFD_Dirichlet_Generate2D_identity) {
    testWGravity(1.1, "identity", "Dirichlet", "Generate2D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_Linear1_MFD_Dirichlet_Generate2D_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_MFD_Dirichlet_Generate2D_diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
#endif
#if NLFV  
  TEST(Analytic00_Linear1_NLFV_Dirichlet_Generate2D_identity) {
    test("identity", "Dirichlet", "Generate2D",
         "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
   TEST(Analytic00_LinearGravity1_NLFV_Dirichlet_Generate2D_identity) {
     testWGravity(1.1, "identity", "Dirichlet", "Generate2D",
                  "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_Linear1_NLFVBFace_Dirichlet_Generate2D_identity) {
     test("identity", "Dirichlet", "Generate2D",
          "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_Dirichlet_Generate2D_identity) {
     testWGravity(1.1, "identity", "Dirichlet", "Generate2D",
                  "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11);
  }
  TEST(Analytic00_Linear1_NLFV_Dirichlet_Generate2D_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_LinearGravity1_NLFV_Dirichlet_Generate2D_diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_Linear1_NLFVBFace_Dirichlet_Generate2D_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-10);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_Dirichlet_Generate2D_diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-10);
  }
#endif


  //  
  // change PC to an assembling PC, change BCs
  //
  // Exact solution
#if ASSEMBLING
  
#if FV  
  TEST(Analytic00_Linear1_FV_Dirichlet_Generate2D_ILU) {
    test("ifpack2: ILUT", "Dirichlet", "Generate2D",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_LinearGravity1_FV_Dirichlet_Generate2D_ILU) {
    testWGravity(1.1, "ifpack2: ILUT", "Dirichlet", "Generate2D",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
#endif
#if MFD  
  TEST(Analytic00_Linear1_MFD_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravity1_MFD_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
#endif
#if NLFV  
  TEST(Analytic00_Linear1_NLFV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
     testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                  "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
#endif
#endif
  
  //  
  // with uniform kr
  //
  // Exact solution
#if FV
  TEST(Analytic00_Linearkr_FV_Dirichlet_Generate2D_Diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "fv", true, AmanziMesh::Entity_kind::FACE, 1e-12);
  }
  TEST(Analytic00_LinearGravitykr_FV_Dirichlet_Generate2D_Diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "fv", true, AmanziMesh::Entity_kind::FACE, 1e-12);
  }
#endif
#if MFD  
  TEST(Analytic00_Linearkr_MFD_Dirichlet_Generate2D_Diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "mixed upwind", true, AmanziMesh::Entity_kind::FACE);
  }
  TEST(Analytic00_LinearGravitykr_MFD_Dirichlet_Generate2D_Diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "mixed upwind", true, AmanziMesh::Entity_kind::FACE);
  }
#endif
#if NLFV  
  TEST(Analytic00_Linearkr_NLFV_Dirichlet_Generate2D_diagonal) {
     test("diagonal", "Dirichlet", "Generate2D",
          "nlfv", true, AmanziMesh::Entity_kind::FACE, 1e-12);
  }
  TEST(Analytic00_LinearGravitykr_NLFV_Dirichlet_Generate2D_diagonal) {
     testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                  "nlfv", true, AmanziMesh::Entity_kind::FACE, 1e-12);
  }
  TEST(Analytic00_Linearkr_NLFVBFace_Dirichlet_Generate2D_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::FACE, 1e-10);
  }
  TEST(Analytic00_LinearGravitykr_NLFVBFace_Dirichlet_Generate2D_diagonal) {
    testWGravity(1.1, "diagonal", "Dirichlet", "Generate2D",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::FACE, 1e-10);
  }
#endif

  //  
  // with uniform but non-1 tensor K
  //
  // Exact solution
#if FV
  TEST(Analytic00_LinearK_FV_DirichletNeumannBox_Generate2D_diagonal) {
    test("diagonal", "DirichletNeumannBox", "Generate2D",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  }
  TEST(Analytic00_LinearGravityK_FV_DirichletNeumannBox_Generate2D_diagonal) {
    testWGravity(1.1, "diagonal", "DirichletNeumannBox", "Generate2D",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  }
#endif
#if MFD  
  TEST(Analytic00_LinearK_MFD_DirichletNeumannBox_Generate2D_diagona) {
    test("diagonal", "DirichletNeumannBox", "Generate2D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
  TEST(Analytic00_LinearGravityK_MFD_DirichletNeumannBox_Generate2D_diagona) {
    testWGravity(1.1, "diagonal", "DirichletNeumannBox", "Generate2D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN);
  }
#endif
#if NLFV  
  TEST(Analytic00_LinearK_NLFV_DirichletNeumannBox_Generate2D_diagona) {
    test("diagonal", "DirichletNeumannBox", "Generate2D",
         "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_LinearGravityK_NLFV_DirichletNeumannBox_Generate2D_diagona) {
    testWGravity(1.1, "diagonal", "DirichletNeumannBox", "Generate2D",
                 "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12);
  }
  TEST(Analytic00_LinearK_NLFVBFace_DirichletNeumannBox_Generate2D_diagona) {
    test("diagonal", "DirichletNeumannBox", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 2e-10);
  }
  TEST(Analytic00_LinearGravityK_NLFVBFace_DirichletNeumannBox_Generate2D_diagona) {
    testWGravity(1.1, "diagonal", "DirichletNeumannBox", "Generate2D",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 2e-10);
  }
#endif

  //  
  // with wiggled mesh
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
#if ASSEMBLING
#if FV
  TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 2.e-2);
  }
  TEST(Analytic00_LinearGravity1_FV_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 2.e-2);
  }
#endif
#if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic00_Linear1_MFD_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 2e-2);
  }
  TEST(Analytic00_LinearGravity1_MFD_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 2e-2);
  }
#endif
#if NLFV
  TEST(Analytic00_Linear1_NLFV_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
         "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "00", 15);
  }
  TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
                 "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "00", 15);
  }
  TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "00", 15);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumannBox_Wiggled_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/random10.exo",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "00", 15);
  }
#endif

  //  
  // with polyhedral mesh
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
#if FV
  TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 3e-2);
  }
  TEST(Analytic00_LinearGravity1_FV_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 3e-2);
  }
#endif
#if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic00_Linear1_MFD_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 2.7e-2);
  }
  TEST(Analytic00_LinearGravity1_MFD_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 2.8e-2);
  }
#endif
#if NLFV
  TEST(Analytic00_Linear1_NLFV_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
         "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 1, "00", 15);
  }
  TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
                 "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "00", 15);
  }
  TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "00", 15);
  }
  TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumannBox_Poly_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/median7x8_filtered.exo",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "00", 15);
  }
#endif

  //
  // Analytic00_Quadratic: test case for quadratic equations
  // polynomial with coefficient=1
  //
#if FV
  TEST(Analytic00_Quadratic1_FV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
  TEST(Analytic00_QuadraticGravity1_FV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
#endif
#if MFD
  TEST(Analytic00_Quadratic1_MFD_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
  TEST(Analytic00_QuadraticGravity1_MFD_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
#endif
#if NLFV
  TEST(Analytic00_Quadratic1_NLFV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
  TEST(Analytic00_QuadraticGravity1_NLFV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
     testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                  "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
  TEST(Analytic00_Quadratic1_NLFVBFace_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
  TEST(Analytic00_QuadraticGravity1_NLFVBFace_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(1.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 0.006, 2);
  }
#endif

  //
  // Analytic02: Tensor diffusion.  Note FV cannot solve this problem, but MFD
  // can (and maybe NLFV as well?)
  //
  // This test replaces old operator_diffusion.cc mixed tests.
#if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic02_MFD_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02");
  }
  TEST(Analytic02_Gravity_MFD_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 2, "02");
  }

  // the original test tested on circle-quad
  TEST(Analytic02_MFD_DirichletNeumannBox_CircleQuad_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/circle_quad10.exo",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02");
  }
  TEST(Analytic02_Gravity_MFD_DirichletNeumannBox_CircleQuad_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/circle_quad10.exo",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02");
  }

  // test in 3D
  TEST(Analytic02_MFD_DirichletNeumannBox_Generate3D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate3D",
         "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 2, "02");
  }
  TEST(Analytic02_Gravity_MFD_DirichletNeumannBox_Generate3D_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate3D",
                 "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02");
  }
#endif

  //
  // With NLFV
  //
#if NLFV
  TEST(Analytic02_NLFV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02", 20);
  }
  TEST(Analytic02_Gravity_NLFV_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02", 20);
  }

  // the original test tested on circle-quad
  TEST(Analytic02_NLFV_DirichletNeumannBox_CircleQuad_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/circle_quad10.exo",
         "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "02", 50);
  }
  TEST(Analytic02_Gravity_NLFV_DirichletNeumannBox_CircleQuad_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/circle_quad10.exo",
                 "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "02", 50);
  }

  // test in 3D
  //
  // These probably should pass?  Not sure why they don't.
  // TEST(Analytic02_NLFV_DirichletNeumannBox_Generate3D_ifpack2_ILUT) {
  //   test("ifpack2: ILUT", "DirichletNeumannBox", "Generate3D",
  //        "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "02", 20);
  // }
  // TEST(Analytic02_Gravity_NLFV_DirichletNeumannBox_Generate3D_ifpack2_ILUT) {
  //   auto ana = Teuchos::rcp(new Analytic02(3, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 0.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate3D",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }

  //
  // With NLFV and Boundary Faces
  // 
  TEST(Analytic02_NLFVwithBndFaces_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02", 15);
  }
  TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumannBox_Generate2D_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate2D",
                 "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-12, 1, "02", 15);
  }

  // the original test tested on circle-quad
  //
  // These probably should pass?  Not sure why the don't.
  TEST(Analytic02_NLFVwithBndFaces_DirichletNeumannBox_CircleQuad_ifpack2_ILUT) {
    test("ifpack2: ILUT", "DirichletNeumannBox", "test/circle_quad10.exo",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-10, 1, "02", 50);
  }
  TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumannBox_CircleQuad_ifpack2_ILUT) {
    testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "test/circle_quad10.exo",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-10, 1, "02", 50);
  }

  // test in 3D
  // TEST(Analytic02_NLFVwithBndFaces_DirichletNeumannBox_Generate3D_ifpack2_ILUT) {
  //   test("ifpack2: ILUT", "DirichletNeumannBox", "GAnalytic02_NLFV_DirichletNeumannBox_Generate3D_ifpack2_ILUTenerate3D",
  //        "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 1, "02", 50);
  // }
  // TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumannBox_Generate3D_ifpack2_ILUT) {
  //   testWGravity(0.1, "ifpack2: ILUT", "DirichletNeumannBox", "Generate3D",
  //                "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1e-11, 1, "02", 50);
  // }
#endif 
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
  TEST(Analytic03b_Linear1_FV_Dirichlet_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "Dirichlet", "Generate2D_HiRes",
         "fv", false, AmanziMesh::Entity_kind::FACE, 2.7e-2, 1, "03");
  }
#endif
#if MFD
  TEST(Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "Dirichlet", "Generate2D_HiRes",
         "mixed upwind", false, AmanziMesh::Entity_kind::FACE, 2.7e-2, 2, "03");
  }
#endif
#if NLFV
  TEST(Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "Dirichlet", "Generate2D",
        "nlfv", false, AmanziMesh::Entity_kind::FACE, 3e-2, 2, "03", 10);
  }
  TEST(Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_ILUT) {
    test("ifpack2: ILUT", "Dirichlet", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::FACE, 3e-2, 2, "03", 10);
  }
#endif

#if FV
  TEST(Analytic03b_Linear1_FV_Dirichlet_Poly_ifpack2_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "fv", false, AmanziMesh::Entity_kind::FACE, 2.7e-2, 1, "03");
  }
#endif
#if MFD
  // on mfd: default, these have a tolerance of 1.e-12.  On TPFA, it is the same as FV?
  TEST(Analytic03b_Linear1_MFD_Dirichlet_Poly_ifpack2_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "mixed upwind", false, AmanziMesh::Entity_kind::FACE, 2.7e-2, 1, "03");
  }
#endif
#if NLFV
  TEST(Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
        "nlfv", false, AmanziMesh::Entity_kind::FACE, 3e-2, 2, "03", 10);
  }
  TEST(Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_diagonal) {
    test("diagonal", "Dirichlet", "Generate2D",
         "nlfv with bfaces", false, AmanziMesh::Entity_kind::FACE, 3e-2, 2, "03", 10);
  }
#endif

}

