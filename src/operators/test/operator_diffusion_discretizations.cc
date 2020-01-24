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
#include "MeshFactory.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
#include "Analytic02.hh"


#include "OperatorDefs.hh"
#include "PDE_DiffusionFV.hh"
// #include "PDE_DiffusionFVwithGravity.hh"
// #include "PDE_DiffusionMFD.hh"
// #include "PDE_DiffusionMFDwithGravity.hh"
// #include "PDE_DiffusionNLFV.hh"
// #include "PDE_DiffusionNLFVwithGravity.hh"
// #include "PDE_DiffusionNLFVwithBndFaces.hh"
// #include "PDE_DiffusionNLFVwithBndFacesGravity.hh"
// #include "UpwindSecondOrder.hh"
// #include "Verification.hh"

#include "DiffusionFixture.hh"

template<class PDE_Diffusion_type>
void test(const Teuchos::RCP<AnalyticBase>& ana,
          const std::string& pc_type,
          const std::string& bc_type,
          const std::string& mesh_type,
          const std::string& disc_type,
          bool symmetric,
          AmanziMesh::Entity_kind scalar_coef, 
          double tol,
          int niters=1) {
  DiffusionFixture fix(ana, mesh_type);
  std::cout << std::endl
            << "Diffusion Test (np=" << fix.comm->getSize() << "): " << disc_type << ", "
            << ana->name() << ", " << pc_type << ", " << mesh_type << std::endl
            << "--------------------------------------------------------------------------------"
            << std::endl;
  fix.discretize<PDE_Diffusion_type, AmanziMesh::FACE>(disc_type);
  std::cout<<"Discretize done"<<std::endl;
  if (scalar_coef != AmanziMesh::Entity_kind::UNKNOWN) fix.scalarCoefficient(scalar_coef);
  if (bc_type == "Dirichlet") {
    fix.setBCsDirichlet();
  } else if (bc_type == "DirichletNeumannBox") {
    fix.setBCsDirichletNeumannBox();
  } else {
    AMANZI_ASSERT(false);
  }      
  fix.setup(pc_type, symmetric);
  for (int i=0; i!=niters-1; ++i) fix.go(0.0);
  fix.go(tol);
}


template<class PDE_Diffusion_type>
void testWGravity(const Teuchos::RCP<AnalyticBase>& ana,
                  double gravity,
                  const std::string& pc_type,
                  const std::string& bc_type,
                  const std::string& mesh_type,
                  const std::string& disc_type,
                  bool symmetric,
                  AmanziMesh::Entity_kind scalar_coef,
                  double tol,
                  int niters=1) {
  DiffusionFixture fix(ana, mesh_type);
  std::cout << std::endl
            << "DiffusionWithGravity Test (np=" << fix.comm->getSize() << "): "
            << disc_type << ", " << ana->name() << ", " << pc_type << ", " << mesh_type << std::endl
            << "--------------------------------------------------------------------------------"
            << std::endl;
  fix.discretizeWithGravity<PDE_Diffusion_type, AmanziMesh::FACE>(disc_type, gravity);
  if (scalar_coef != AmanziMesh::Entity_kind::UNKNOWN) fix.scalarCoefficient(scalar_coef);
  if (bc_type == "Dirichlet") {
    fix.setBCsDirichlet();
  } else if (bc_type == "DirichletNeumannBox") {
    fix.setBCsDirichletNeumannBox();
  } else {
    AMANZI_ASSERT(false);
  }      
  fix.setup(pc_type, symmetric);
  for (int i=0; i!=niters-1; ++i) fix.go(0.0);
  fix.go(tol);
}


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
  //    preconditioner (diagonal, HYPRE)
  //    analytic problem:
  //            Analytic00_Linear: coef=1, linear solution (exact for all methods!)
  //            Analytic00_Quadratic: coef=1, quadratic solution (exact for second order methods!)

  //            Analytic00b: 3D variant of Analytic00
  //            Analytic03: discontinuous scalar Tensor coef, no scalar coef
  //            Analytic03b: same as 03, but with Tensor coef=1 and scalar coef as the scalar
  //
  // Not all combinations make sense, obviously.  2D problems must be run with
  // 2D meshes, and 3D with 3D.

  TEST(ConstructorDestructor) {
    auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
    DiffusionFixture fix(ana, "Generate2D");
    fix.discretize<Operators::PDE_DiffusionFV, AmanziMesh::FACE>("fv");
    fix.setBCsDirichlet();
    fix.setup("diagonal", true);

    fix.global_op->Zero();
    fix.op->UpdateMatrices(Teuchos::null, fix.solution.ptr());

    CompositeVector& rhs = *fix.global_op->rhs();
    {
      auto rhs_c = rhs.ViewComponent<AmanziDefaultHost>("cell", false);
      for (int c=0; c!=fix.mesh->num_entities(AmanziMesh::Entity_kind::CELL,
              AmanziMesh::Parallel_type::OWNED); ++c) {
        const auto& xc = fix.mesh->cell_centroid(c);
        rhs_c(c,0) += ana->source_exact(xc, 0.0) * fix.mesh->cell_volume(c);
      }
    }

    fix.op->ApplyBCs(true, true, true);
    fix.global_op->UpdatePreconditioner();
    

  }

  
  //
  // Analytic00_Linear: tests exactness of no-gravity case for linear
  // polynomial with coefficient=1
  //
  // Exact solution
  TEST(Analytic00_Linear1_FV_Dirichlet_Generate2D_diagonal) {
    auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
    test<Operators::PDE_DiffusionFV>(
        ana, "diagonal", "Dirichlet", "Generate2D",
        "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  }

  // TEST(Analytic00_LinearGravity1_FV_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "diagonal", "Dirichlet", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_MFD_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "diagonal", "Dirichlet", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_MFD_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "diagonal", "Dirichlet", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_NLFV_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "diagonal", "Dirichlet", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_NLFV_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "diagonal", "Dirichlet", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_NLFVBFace_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "diagonal", "Dirichlet", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_NLFVBFace_Dirichlet_Generate2D_diagonal) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "diagonal", "Dirichlet", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }

  //  
  // change PC, change BCs
  //
  // Exact solution
  // TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }


  //  
  // with uniform kr
  //
  // Exact solution
  // TEST(Analytic00_Linearkr_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionFV>(
  //       ana, "HypreAMG", "Dirichlet", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::FACE, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravitykr_FV_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "Dirichlet", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::FACE, 1.e-12);
  // }
  // TEST(Analytic00_Linearkr_MFD_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "Dirichlet", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::CELL, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravitykr_MFD_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "HypreAMG", "Dirichlet", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::CELL, 1.e-12);
  // }
  // TEST(Analytic00_Linearkr_NLFV_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "Dirichlet", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::FACE, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravitykr_NLFV_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "Dirichlet", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::FACE, 1.e-12);
  // }
  // TEST(Analytic00_Linearkr_NLFVBFace_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "Dirichlet", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::FACE, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravitykr_NLFVBFace_Dirichlet_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 3.14, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "HypreAMG", "Dirichlet", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::FACE, 1.e-12);
  // }


  //  
  // with uniform K
  //
  // Exact solution
  // TEST(Analytic00_LinearK_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 0.0));
  //   test<Operators::PDE_DiffusionFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravityK_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearK_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravityK_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearK_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravityK_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearK_NLFVBFace_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravityK_NLFVBFace_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 3.14, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }


  //  
  // with wiggled mesh
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
  // TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 2.e-2);
  // }
  // TEST(Analytic00_LinearGravity1_FV_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 2.e-2);
  // }
  // TEST(Analytic00_Linear1_MFD_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_MFD_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_NLFV_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumannBox_Wiggled_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/random10.exo",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }


  //  
  // with polyhedral mesh
  //
  // NOTE: FV cannot be exact
  //       MFD is still exact
  //       NLFV must converge to be exact.  Also, NLFV loses symmetry as it iterates.
  // TEST(Analytic00_Linear1_FV_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 2.7e-2);
  // }
  // TEST(Analytic00_LinearGravity1_FV_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, 2.8e-2);
  // }
  // TEST(Analytic00_Linear1_MFD_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_LinearGravity1_MFD_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic00_Linear1_NLFV_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic00_LinearGravity1_NLFV_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic00_LinearGravity1_NLFVBFace_DirichletNeumannBox_Poly_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(1, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "test/median7x8_filtered.exo",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }


  //
  // Analytic00_Quadratic: test case for quadratic equations
  // polynomial with coefficient=1
  //
  // TEST(Analytic00_Quadratic1_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_QuadraticGravity1_FV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "fv", true, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_Quadratic1_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_QuadraticGravity1_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_Quadratic1_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_QuadraticGravity1_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", true, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_Quadratic1_NLFVBFace_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 0.0));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }
  // TEST(Analytic00_QuadraticGravity1_NLFVBFace_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic00(2, 1.0, 1.0, 1.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 1.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, .006);
  // }

  //
  // Analytic02: Tensor diffusion.  Note FV cannot solve this problem, but MFD
  // can (and maybe NLFV as well?)
  //
  // This test replaces old operator_diffusion.cc mixed tests.
  // TEST(Analytic02_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic02_Gravity_MFD_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2, 0.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }

  // the original test tested on circle-quad
  // TEST(Analytic02_MFD_DirichletNeumannBox_CircleQuad_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/circle_quad10.exo",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic02_Gravity_MFD_DirichletNeumannBox_CircleQuad_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2, 0.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "test/circle_quad10.exo",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }

  // test in 3D
  // TEST(Analytic02_MFD_DirichletNeumannBox_Generate3D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(3));
  //   test<Operators::PDE_DiffusionMFD>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate3D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }
  // TEST(Analytic02_Gravity_MFD_DirichletNeumannBox_Generate3D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(3, 0.1));
  //   testWGravity<Operators::PDE_DiffusionMFDwithGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "Generate3D",
  //       "mixed", true, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12);
  // }


  //
  // With NLFV
  // 
  // TEST(Analytic02_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic02_Gravity_NLFV_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }

  // the original test tested on circle-quad
  // TEST(Analytic02_NLFV_DirichletNeumannBox_CircleQuad_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/circle_quad10.exo",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }
  // TEST(Analytic02_Gravity_NLFV_DirichletNeumannBox_CircleQuad_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "test/circle_quad10.exo",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }

  // test in 3D
  //
  // These probably should pass?  Not sure why they don't.
  // TEST(Analytic02_NLFV_DirichletNeumannBox_Generate3D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(3));
  //   test<Operators::PDE_DiffusionNLFV>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate3D",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }
  // TEST(Analytic02_Gravity_NLFV_DirichletNeumannBox_Generate3D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(3, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "Generate3D",
  //       "nlfv", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }

  //
  // With NLFV and Boundary Faces
  // 
  // TEST(Analytic02_NLFVwithBndFaces_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }
  // TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumannBox_Generate2D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "Generate2D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-12, 20);
  // }

  // the original test tested on circle-quad
  //
  // These probably should pass?  Not sure why the don't.
  // TEST(Analytic02_NLFVwithBndFaces_DirichletNeumannBox_CircleQuad_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "test/circle_quad10.exo",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }
  // TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumannBox_CircleQuad_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(2, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "test/circle_quad10.exo",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }

  // test in 3D
  // TEST(Analytic02_NLFVwithBndFaces_DirichletNeumannBox_Generate3D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(3));
  //   test<Operators::PDE_DiffusionNLFVwithBndFaces>(
  //       ana, "HypreAMG", "DirichletNeumannBox", "Generate3D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }
  // TEST(Analytic02_Gravity_NLFVwithBndFaces_DirichletNeumannBox_Generate3D_HypreAMG) {
  //   auto ana = Teuchos::rcp(new Analytic02(3, 0.1));
  //   testWGravity<Operators::PDE_DiffusionNLFVwithBndFacesGravity>(
  //       ana, 0.1, "HypreAMG", "DirichletNeumannBox", "Generate3D",
  //       "nlfv with bfaces", false, AmanziMesh::Entity_kind::UNKNOWN, 1.e-11, 50);
  // }
  
}
