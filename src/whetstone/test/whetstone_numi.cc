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
#include "Point.hh"

#include "NumericalIntegration.hh"
#include "Polynomial.hh"

class TestFunction : public Amanzi::WhetStone::WhetStoneFunction {
 public:
  virtual double Value(const Amanzi::AmanziGeometry::Point& xp) const
  {
    double a = norm(xp - Amanzi::AmanziGeometry::Point(0.5, 0.5));
    return (a < 0.25) ? 1.0 : 0.0;
    // return (1.0 + tanh((a - 0.5) * 10)) / 2;
  }
};


/* **************************************************************** */
TEST(NUMI_CELL_2D_EULER_FORMULA)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: Euler's formula for polygon" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList()); 
  fac_list->set<bool>("request edges", true); 
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  NumericalIntegration numi(mesh);

  int cell(0);
  double val;

  // 0th-order polynomial
  Polynomial poly(2, 0);
  poly(0, 0) = 1.0;
  val = numi.IntegratePolynomialCell(cell, poly);

  printf("order=0  value=%10.6g\n", val);
  CHECK_CLOSE(val, mesh->getCellVolume(cell), 1e-10);

  // 1st-order polynomial
  poly.Reshape(2, 1);
  poly(1, 0) = 2.0;
  poly(1, 1) = 3.0;
  poly.set_origin(mesh->getCellCentroid(cell));
  val = numi.IntegratePolynomialCell(cell, poly);

  printf("order=1  value=%10.6g\n", val);
  CHECK_CLOSE(val, mesh->getCellVolume(cell), 1e-10);
}


/* **************************************************************** */
TEST(NUMI_CELL_2D_QUADRATURE_POLYGON)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: quadrature rules in 2D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList()); 
  fac_list->set<bool>("request edges", true); 
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  NumericalIntegration numi(mesh);

  int cell(0), face(0);
  double val1, val2;

  // linear polynomial centered at the origin
  Polynomial poly(2, 1);
  poly(1, 0) = 2.0;
  poly(1, 1) = 3.0;
  poly.set_origin(mesh->getCellCentroid(cell));

  std::vector<const WhetStoneFunction*> polys(1);
  polys[0] = &poly;

  for (int order = 1; order < 14; ++order) {
    val1 = numi.IntegrateFunctionsTriangulatedCell(cell, polys, order);

    printf("order=%d  value=%10.6g\n", order, val1);
    CHECK_CLOSE(val1, poly.Value(mesh->getCellCentroid(cell)), 1e-12);
  }

  // cross-comparison of integrators
  for (int order = 1; order < 14; ++order) {
    poly.Reshape(2, order, true);
    poly(order, 0) = 1.0;

    val1 = numi.IntegrateFunctionsTriangulatedCell(cell, polys, order);
    val2 = numi.IntegratePolynomialCell(cell, poly);

    printf("CELL: order=%d  values: %10.6g %10.6g\n", order, val1, val2);
    CHECK_CLOSE(val1, val2, 1e-12);

    // face of a polygon
    val1 = numi.IntegrateFunctionsTriangulatedFace(face, polys, order);
    val2 = numi.IntegratePolynomialFace(face, poly);
    printf("FACE: order=%d  values: %10.6g %10.6g\n", order, val1, val2);
    CHECK_CLOSE(val1, val2, 1e-12);
  }
}


/* **************************************************************** */
TEST(NUMI_CELL_2D_QUADRATURE_SQUARE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: quadrature rule for square" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 2, 2);

  NumericalIntegration numi(mesh);

  int cell(3), face(9);
  double val, exact;

  for (int order = 1; order < 14; ++order) {
    Polynomial poly(2, order);
    poly(order, 0) = 1.0;

    std::vector<const WhetStoneFunction*> polys(1);
    polys[0] = &poly;

    val = numi.IntegrateFunctionsTriangulatedCell(cell, polys, poly.order());
    exact = 1.0 / (order + 1);

    printf("CELL: order=%d  value=%10.6g  exact=%10.6g\n", order, val, exact);
    CHECK_CLOSE(val, exact, 1e-12);

    // face of a square
    val = numi.IntegrateFunctionsTriangulatedFace(face, polys, order);
    printf("FACE: order=%d  values: %10.6g %10.6g  err=%10.6g\n", order, val, exact, val - exact);
    CHECK_CLOSE(val, exact, 1e-12);

    // add y-component to polynomial
    int k = MonomialSpaceDimension(2, order);
    poly(order, k - 1) = 1.0;
    val = numi.IntegrateFunctionsTriangulatedCell(cell, polys, order);
    CHECK_CLOSE(val, 2 * exact, 1e-12);

    // functions with sharp features
    // -- quadrature based on partition
    TestFunction f;
    if (order > 6) {
      polys[0] = &f;
      val = numi.IntegrateFunctionsTriangulatedCell(cell, polys, order);
      CHECK_CLOSE(val, 0.1963, 0.1);
    }

    // -- quadrature based on tensor product
    int mmax = order / 2;
    double val2(0.0);
    exact = 0.196349540849362;
    AmanziGeometry::Point q(2);

    for (int i = 0; i <= mmax; ++i) {
      q[0] = q1d_points[mmax][i];
      for (int j = 0; j <= mmax; ++j) {
        q[1] = q1d_points[mmax][j];
        val2 += f.Value(q) * q1d_weights[mmax][i] * q1d_weights[mmax][j];
      }
    }
    std::cout << "n=" << order << " " << val << " " << val2 << " errs: " << val - exact << " "
              << val2 - exact << std::endl;
  }
}


/* **************************************************************** */
TEST(NUMI_CELL_3D_QUADRATURE_POLYHEDRON)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: quadrature rules in 3D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList()); 
  fac_list->set<bool>("request edges", true); 
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/dodecahedron.exo");

  NumericalIntegration numi(mesh);

  int cell(0), face(3);
  double val1, val2;

  // linear polynomial centered at the origin
  Polynomial poly(3, 1);
  poly(1, 0) = 2.0;
  poly(1, 1) = 3.0;
  poly(1, 2) = 4.0;
  poly.set_origin(mesh->getCellCentroid(cell));

  std::vector<const WhetStoneFunction*> polys(1);
  polys[0] = &poly;

  for (int order = 1; order < 7; ++order) {
    val1 = numi.IntegrateFunctionsTriangulatedCell(cell, polys, order);

    printf("order=%d  value=%10.6g\n", order, val1);
    CHECK_CLOSE(val1, poly.Value(mesh->getCellCentroid(cell)), 1e-12);
  }

  // cross-comparison of integrators
  for (int order = 1; order < 7; ++order) {
    poly.Reshape(3, order, true);
    poly(order, 0) = 1.0;

    val1 = numi.IntegrateFunctionsTriangulatedCell(cell, polys, order);
    val2 = numi.IntegratePolynomialCell(cell, poly);

    printf("CELL: order=%d  value= %12.6f  diff=%10.6g\n", order, val1, val1 - val2);
    CHECK_CLOSE(val1, val2, 5e-11 * std::max(1.0, std::fabs(val1)));

    // face of a polyhedron
    val1 = numi.IntegrateFunctionsTriangulatedFace(face, polys, order);
    val2 = numi.IntegratePolynomialFace(face, poly);
    printf("FACE: order=%d  values: %10.6g %10.6g\n", order, val1, val2);
    CHECK_CLOSE(val1, val2, 1e-12 * std::max(1.0, std::fabs(val1)));
  }
}


/* ******************************************************************
* 3D cube: Amanzi mesh
****************************************************************** */
TEST(NUMI_CELL_3D_QUADRATURE_CUBE)
{
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: quadrature rule for cube" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  AmanziMesh::MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(AmanziMesh::Preference({ AmanziMesh::Framework::MSTK }));
  Teuchos::RCP<AmanziMesh::Mesh> mesh =
    meshfactory.create(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 2, 2, 2);

  NumericalIntegration numi(mesh);

  int cell(7), face(31);
  double val, exact;

  for (int order = 1; order < 7; ++order) {
    Polynomial poly(3, order);
    poly(order, 0) = 1.0;

    std::vector<const WhetStoneFunction*> polys(1);
    polys[0] = &poly;

    val = numi.IntegrateFunctionsTriangulatedCell(cell, polys, poly.order());
    exact = 1.0 / (order + 1);

    printf("CELL: order=%d  value=%10.6g  exact=%10.6g\n", order, val, exact);
    CHECK_CLOSE(val, exact, 1e-13);

    // face of a polyhedron
    val = numi.IntegrateFunctionsTriangulatedFace(face, polys, poly.order());
    printf("FACE: order=%d  values: %10.6g %10.6g\n", order, val, exact);
    CHECK_CLOSE(val, exact, 1e-12 * std::fabs(val));
  }
}
