/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "GrammMatrix.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_CrouzeixRaviartAnyOrder.hh"
#include "MFD3D_CrouzeixRaviartSerendipity.hh"
#include "MFD3D_Lagrange.hh"
#include "MFD3D_LagrangeAnyOrder.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "Tensor.hh"


/* ******************************************************************
* Crouzier-Raviart 2D and 3D elements
****************************************************************** */
void
HighOrderCrouzeixRaviart(int dim, std::string file_name)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: High-order Crouzeix Raviart in " << dim << "D"
            << " file=" << file_name << std::endl;
  auto comm = Amanzi::getDefaultComm();
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
  if (dim == 3) fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, Teuchos::null, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(file_name);

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);
  MFD3D_CrouzeixRaviart mfd_lo(plist, mesh);
  MFD3D_CrouzeixRaviartAnyOrder mfd_ho(plist, mesh);

  int cell(0);
  DenseMatrix G1, A1, Ak;

  Tensor T(dim, 2);
  T(0, 0) = T(1, 1) = 2.0;
  T(0, 1) = T(1, 0) = 1.0;
  T(dim - 1, dim - 1) = 2.0;

  // 1st-order scheme
  mfd_lo.set_order(1);
  mfd_lo.StiffnessMatrix(cell, T, A1);

  // 1st-order scheme (new algorithm)
  mfd_ho.set_order(1);
  mfd_ho.StiffnessMatrix(cell, T, Ak);

  printf("Stiffness (sub)matrix for order = 1\n");
  PrintMatrix(A1, "%8.4f ");

  A1 -= Ak;
  CHECK(A1.NormInf() <= 1e-10);

  // high-order scheme (new algorithm)
  int kmax = (dim == 3) ? 3 : 4;
  for (int k = 2; k < kmax; ++k) {
    mfd_ho.set_order(k);
    mfd_ho.StiffnessMatrix(cell, T, Ak);

    printf("Stiffness (sub)matrix for order=%d, size=%d\n", k, Ak.NumRows());
    PrintMatrix(Ak, "%8.4f ", 12);

    // verify SPD propery
    int nrows = Ak.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(Ak(i, i) > 0.0);

    // verify exact integration property
    const DenseMatrix& G = mfd_ho.G();

    PolynomialOnMesh integrals;
    NumericalIntegration numi(mesh);
    numi.UpdateMonomialIntegralsCell(cell, 2 * k, integrals);

    Polynomial ptmp, poly(dim, k);
    Basis_Regularized basis;
    basis.Init(mesh, cell, k, ptmp);

    GrammMatrixGradients(T, poly, integrals, basis, G1);

    G1(0, 0) = 1.0;
    G1.Inverse();
    G1(0, 0) = 0.0;
    G1 -= G;
    CHECK(G1.NormInf() <= 1e-10);
  }
}


TEST(HIGH_ORDER_CROUZEIX_RAVIART)
{
  HighOrderCrouzeixRaviart(2, "test/one_pentagon.exo");
  HighOrderCrouzeixRaviart(3, "test/cube_unit.exo");
}


/* ******************************************************************
* Incorrect Serendipity Crouzier-Raviart 2D element (for testing)
****************************************************************** */
void
HighOrderCrouzeixRaviartSerendipity(int dim, std::string file_name)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order Crouzeix Raviart Serendipity element in " << dim
            << "D, file=" << file_name << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
  if (dim == 3) fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(file_name);

  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);
  MFD3D_CrouzeixRaviartSerendipity mfd(plist, mesh);

  for (int c = 0; c < ncells; ++c) {
    if (mesh->getCellNumFaces(c) < 4) continue;

    Tensor T(dim, 1);
    T(0, 0) = 1.0;

    // high-order schemes
    int kmax = (dim == 3) ? 3 : 4;
    for (int k = 1; k < kmax; ++k) {
      DenseMatrix G1, Ak;

      mfd.set_order(k);
      mfd.StiffnessMatrix(c, T, Ak);

      printf("Stiffness (sub)matrix for order=%d, cell=%d, size=%d\n", k, c, Ak.NumRows());
      PrintMatrix(Ak, "%8.3f ", 12);

      // verify SPD propery
      int nrows = Ak.NumRows();
      for (int i = 0; i < nrows; i++) CHECK(Ak(i, i) > 0.0);

      // verify exact integration property
      const DenseMatrix& G = mfd.G();

      PolynomialOnMesh integrals;
      NumericalIntegration numi(mesh);
      numi.UpdateMonomialIntegralsCell(c, 2 * k, integrals);

      Polynomial ptmp, poly(dim, k);
      Basis_Regularized basis;
      basis.Init(mesh, c, k, ptmp);

      GrammMatrixGradients(T, poly, integrals, basis, G1);

      G1(0, 0) = 1.0;
      G1.Inverse();
      G1(0, 0) = 0.0;
      G1 -= G;
      std::cout << "Norm of dG=" << G1.NormInf() << " " << G.NormInf() << std::endl;
      CHECK(G1.NormInf() <= 1e-11 * G.NormInf());
    }
  }
}


TEST(HIGH_ORDER_CROUZEIX_RAVIART_SERENDIPITY)
{
  HighOrderCrouzeixRaviartSerendipity(2, "test/two_cell2_dist.exo");
  HighOrderCrouzeixRaviartSerendipity(2, "test/one_pentagon.exo");
  HighOrderCrouzeixRaviartSerendipity(3, "test/cube_unit.exo");
}


/* ******************************************************************
* Lagrange 2D element
****************************************************************** */
void
HighOrderLagrange2D(std::string file_name)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order Lagrange element, file=" << file_name << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
  fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(file_name);

  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);

  MFD3D_Lagrange mfd_lo(plist, mesh);
  MFD3D_LagrangeAnyOrder mfd_ho(plist, mesh);

  for (int c = 0; c < ncells; ++c) {
    DenseMatrix G1, A1, Ak;

    Tensor T(2, 2);
    T(0, 0) = T(1, 1) = 2.0;
    T(0, 1) = T(1, 0) = 0.1;

    // 1st-order scheme
    mfd_lo.StiffnessMatrix(c, T, A1);

    // 1st-order scheme (new algorithm)
    mfd_ho.set_order(1);
    mfd_ho.StiffnessMatrix(c, T, Ak);

    printf("Stiffness (sub)matrix for order=1, cell=%d\n", c);
    PrintMatrix(Ak, "%8.4f ");

    A1 -= Ak;
    CHECK(A1.NormInf() <= 1e-10);

    // high-order scheme (new algorithm)
    for (int k = 2; k < 5; ++k) {
      mfd_ho.set_order(k);
      mfd_ho.StiffnessMatrix(c, T, Ak);

      printf("Stiffness (sub)matrix for order=%d, cell=%d, size=%d\n", k, c, Ak.NumRows());
      PrintMatrix(Ak, "%8.4f ", 12);

      // verify SPD propery
      int nrows = Ak.NumRows();
      for (int i = 0; i < nrows; i++) CHECK(Ak(i, i) > 0.0);

      // verify exact integration property
      const DenseMatrix& G = mfd_ho.G();

      PolynomialOnMesh integrals;
      NumericalIntegration numi(mesh);
      numi.UpdateMonomialIntegralsCell(c, 2 * k, integrals);

      Polynomial ptmp, poly(mesh->getSpaceDimension(), k);
      Basis_Regularized basis;
      basis.Init(mesh, c, k, ptmp);

      GrammMatrixGradients(T, poly, integrals, basis, G1);

      G1(0, 0) = 1.0;
      G1.Inverse();
      G1(0, 0) = 0.0;
      G1 -= G;
      CHECK(G1.NormInf() <= 1e-12 * G.NormInf());
    }
  }
}


TEST(HIGH_ORDER_LAGRANGE_2D)
{
  HighOrderLagrange2D("test/one_pentagon.exo");
  HighOrderLagrange2D("test/two_cell2_dist.exo");
}


/* ******************************************************************
* Lagrange 3D element
****************************************************************** */
void
HighOrderLagrange3D(const std::string& filename1, const std::string& filename2)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order 3D Lagrange element, comparing: " << filename1 << " and "
            << filename2 << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
  fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh1 = meshfactory.create(filename1);
  Teuchos::RCP<Mesh> mesh2 = meshfactory.create(filename2);

  DenseMatrix G1, A1, A2;
  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);

  Tensor T(3, 1);
  T(0, 0) = 2.0;

  {
    MFD3D_Lagrange mfd_lo(plist, mesh2);
    MFD3D_LagrangeAnyOrder mfd_ho(plist, mesh2);

    // 1st-order scheme
    mfd_lo.StiffnessMatrix(0, T, A1);

    // 1st-order scheme (new algorithm)
    mfd_ho.set_order(1);
    mfd_ho.StiffnessMatrix(0, T, A2);

    printf("Stiffness (sub)matrix for order=1\n");
    PrintMatrix(A2, "%8.4f ");

    A1 -= A2;
    CHECK(A1.NormInf() <= 1e-10 * A2.NormInf());
  }

  // high-order scheme (new algorithm)
  MFD3D_LagrangeAnyOrder mfd1(plist, mesh1);
  MFD3D_LagrangeAnyOrder mfd2(plist, mesh2);
  for (int k = 1; k < 3; ++k) {
    mfd1.set_order(k);
    mfd1.StiffnessMatrix(0, T, A1);

    mfd2.set_order(k);
    mfd2.StiffnessMatrix(0, T, A2);

    printf("Stiffness (sub)matrix for order=%d, size=%d\n", k, A1.NumRows());
    PrintMatrix(A1, "%8.4f ", 12);

    // verify SPD propery
    int nrows = A1.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(A1(i, i) > 0.0);

    A1 -= A2;
    CHECK(A1.NormInf() <= 1e-10 * A2.NormInf());

    // verify exact integration property
    const DenseMatrix& G = mfd2.G();

    PolynomialOnMesh integrals;
    NumericalIntegration numi(mesh2);
    numi.UpdateMonomialIntegralsCell(0, 2 * k, integrals);

    Polynomial ptmp, poly(3, k);
    Basis_Regularized basis;
    basis.Init(mesh2, 0, k, ptmp);

    GrammMatrixGradients(T, poly, integrals, basis, G1);

    G1(0, 0) = 1.0;
    G1.Inverse();
    G1(0, 0) = 0.0;
    G1 -= G;
    CHECK(G1.NormInf() <= 1e-12 * G.NormInf());
  }
}


TEST(HIGH_ORDER_LAGRANGE_3D)
{
  HighOrderLagrange3D("test/cube_unit.exo", "test/cube_unit_rotated.exo");
  HighOrderLagrange3D("test/cube_half.exo", "test/cube_half.exo");
  HighOrderLagrange3D("test/parallepiped.exo", "test/parallepiped_rotated.exo");
}


/* ******************************************************************
* Serendipity 2D and 3D Lagrange elements
****************************************************************** */
void
HighOrderLagrangeSerendipity(const std::string& filename)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order Lagrange Serendipity element, file=" << filename << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
  fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(filename);

  int d = mesh->getSpaceDimension();
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);

  MFD3D_Lagrange mfd_lo(plist, mesh);
  MFD3D_LagrangeSerendipity mfd_ho(plist, mesh);

  for (int c = 0; c < ncells; ++c) {
    if (mesh->getCellNumFaces(c) < 4) continue;

    int rank = (d == 3) ? 1 : 2;
    Tensor T(d, rank);
    T(0, 0) = 2.0;
    if (rank == 2) {
      T(1, 1) = 2.0;
      T(0, 1) = T(1, 0) = 1.0;
    }

    // 1st-order scheme
    DenseMatrix G1, A1, Ak;
    mfd_lo.StiffnessMatrix(c, T, A1);

    mfd_ho.set_order(1);
    mfd_ho.StiffnessMatrix(c, T, Ak);

    printf("Stiffness (sub)matrix for order=1, cell=%d\n", c);
    PrintMatrix(Ak, "%8.4f ");

    A1 -= Ak;
    CHECK(A1.NormInf() <= 1e-10);

    // high-order schemes
    for (int k = 2; k < 4; ++k) {
      mfd_ho.set_order(k);
      mfd_ho.StiffnessMatrix(c, T, Ak);

      printf("Stiffness (sub)matrix for order=%d, cell=%d, size=%d\n", k, c, Ak.NumRows());
      PrintMatrix(Ak, "%8.3f ", 12);

      // verify SPD propery
      int nrows = Ak.NumRows();
      for (int i = 0; i < nrows; i++) CHECK(Ak(i, i) > 0.0);

      // verify exact integration property
      const DenseMatrix& G = mfd_ho.G();

      PolynomialOnMesh integrals;
      NumericalIntegration numi(mesh);
      numi.UpdateMonomialIntegralsCell(0, 2 * k, integrals);

      Polynomial ptmp, poly(d, k);
      Basis_Regularized basis;
      basis.Init(mesh, 0, k, ptmp);

      GrammMatrixGradients(T, poly, integrals, basis, G1);

      G1(0, 0) = 1.0;
      G1.Inverse();
      G1(0, 0) = 0.0;
      G1 -= G;
      std::cout << "Norm of dG=" << G1.NormInf() << " " << G.NormInf() << std::endl;
      CHECK(G1.NormInf() <= 1e-12 * G.NormInf());
    }
  }
}


TEST(HIGH_ORDER_LAGRANGE_SERENDIPITY)
{
  HighOrderLagrangeSerendipity("test/two_cell2_dist.exo");
  HighOrderLagrangeSerendipity("test/one_pentagon.exo");
  HighOrderLagrangeSerendipity("test/cube_unit.exo");
}


/* ******************************************************************
* Surface Lagrange element
****************************************************************** */
TEST(HIGH_ORDER_LAGRANGE_SURFACE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order Lagrange element on surface" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
  MeshFactory meshfactory1(comm, gm, fac_list);
  meshfactory1.set_preference(Preference({ Framework::MSTK }));
  fac_list->set<bool>("request edges", true);
  fac_list->set<bool>("request faces", true);
  MeshFactory meshfactory2(comm, gm, fac_list);
  meshfactory2.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh2d = meshfactory1.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  Teuchos::RCP<Mesh> mesh3d = meshfactory2.create("test/cube_unit.exo");

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 2);

  MFD3D_LagrangeAnyOrder mfd_2d(plist, mesh2d);
  MFD3D_LagrangeAnyOrder mfd_ho(plist, mesh3d);

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  // 1st-order scheme
  DenseMatrix A2d, A3d;
  mfd_ho.StiffnessMatrixSurface(0, T, A3d);

  printf("Stiffness (sub)matrix for order=2, size=%d\n", A3d.NumRows());
  PrintMatrix(A3d, "%9.3f ", 12);

  // compare with flat construction
  mfd_2d.StiffnessMatrix(0, T, A2d);

  A3d -= A2d;
  CHECK(A3d.NormInf() <= 1e-12 * A2d.NormInf());
}
