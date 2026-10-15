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

#define FV 0
#define SO 
#define MFD 0
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



  // TEST(Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_diagonal)
  // {
  //   std::cout << "Test: "
  //             << "Analytic03b_Linear1_NLFV_Dirichlet_Poly_ifpack2_diagonal" << std::endl;
  //   test("diagonal",
  //        "Dirichlet",
  //        "structured2d",
  //        2,
  //        10,
  //        "nlfv",
  //        false,
  //        AmanziMesh::Entity_kind::FACE,
  //        3e-2,
  //        2,
  //        "03",
  //        10);
  // }
  // TEST(Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_diagonal)
  // {
  //   std::cout << "Test: "
  //             << "Analytic03b_Linear1_NLFVBFace_Dirichlet_Poly_ifpack2_diagonal" << std::endl;
  //   test("diagonal",
  //        "Dirichlet",
  //        "structured2d",
  //        2,
  //        10,
  //        "nlfv with bfaces",
  //        false,
  //        AmanziMesh::Entity_kind::FACE,
  //        3e-2,
  //        2,
  //        "03",
  //        10);
  // }

  // TEST(Analytic00_Linear1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT)
  // {
  //   std::cout << "Test: "
  //             << "Analytic00_Linear1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
  //   test("ifpack2: ILUT",
  //        "DirichletNeumann",
  //        "structured2d",
  //        2,
  //        5,
  //        "nlfv",
  //        true,
  //        AmanziMesh::Entity_kind::UNKNOWN,
  //        1e-12);
  // }

  // TEST(Analytic00_Linear1_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT)
  // {
  //   std::cout << "Test: "
  //             << "Analytic00_Linear1_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
  //   test("ifpack2: ILUT",
  //        "DirichletNeumann",
  //        "structured2d",
  //        2,
  //        5,
  //        "nlfv with bfaces",
  //        false,
  //        AmanziMesh::Entity_kind::UNKNOWN,
  //        1e-12);
  // }

  TEST(Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT0)
  {
    std::cout << "Test: "
              << "Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         5,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.004,
         1,
         "02");
  }

  TEST(Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT2)
  {
    std::cout << "Test: "
              << "Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         20,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.004,
         1,
         "02");
  }

    TEST(Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT3)
  {
    std::cout << "Test: "
              << "Analytic02_NLFV_DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
    test("ifpack2: ILUT",
         "Dirichlet",
         "structured2d",
         2,
         40,
         "nlfv",
         false,
         AmanziMesh::Entity_kind::UNKNOWN,
         0.004,
         1,
         "02");
  }

  // TEST(Analytic02_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT0)
  // {
  //   std::cout << "Test: "
  //             << "Analytic02_NLFVwithBnd__DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
  //   test("ifpack2: ILUT",
  //        "DirichletNeumann",
  //        "structured2d",
  //        2,
  //        5,
  //        "nlfv with bfaces",
  //        false,
  //        AmanziMesh::Entity_kind::UNKNOWN,
  //        0.004,
  //        1,
  //        "02");
  // }

  //   TEST(Analytic02_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT1)
  // {
  //   std::cout << "Test: "
  //             << "Analytic02_NLFVwithBnd__DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
  //   test("ifpack2: ILUT",
  //        "DirichletNeumann",
  //        "structured2d",
  //        2,
  //        10,
  //        "nlfv with bfaces",
  //        false,
  //        AmanziMesh::Entity_kind::UNKNOWN,
  //        0.004,
  //        1,
  //        "02");
  // }

  //   TEST(Analytic02_NLFVBFace_DirichletNeumann_structured2d_ifpack2_ILUT2)
  // {
  //   std::cout << "Test: "
  //             << "Analytic02_NLFVwithBnd__DirichletNeumann_structured2d_ifpack2_ILUT" << std::endl;
  //   test("ifpack2: ILUT",
  //        "DirichletNeumann",
  //        "structured2d",
  //        2,
  //        20,
  //        "nlfv with bfaces",
  //        false,
  //        AmanziMesh::Entity_kind::UNKNOWN,
  //        0.004,
  //        1,
  //        "02");
  // }

}
