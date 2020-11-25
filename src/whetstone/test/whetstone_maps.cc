/*
  WhetStone

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

// TPLs
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Mesh.hh"
#include "MeshFactory.hh"

// Amanzi::WhetStone
#include "DenseMatrix.hh"
#include "DG_Modal.hh"
#include "MeshMaps_VEM.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Lagrange.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"


/* ****************************************************************
* Test of determinant of transformation
**************************************************************** */
TEST(DG_MAP_DETERMINANT_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Convergence of determinant." << std::endl;
  auto comm = Amanzi::getDefaultComm();

  // create two meshes
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<Mesh> mesh0 = meshfactory.create("test/one_pentagon.exo");
  Teuchos::RCP<Mesh> mesh1 = meshfactory.create("test/one_pentagon.exo");

  // deform the second mesh
  int dim(2), cell(0), nnodes(5), nfaces(5);
  AmanziGeometry::Point xv(dim);
  Entity_ID_List nodeids, faces;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    xv[0] = (xv[0] + xv[0] * xv[0] + xv[0] * xv[0] * xv[0]) / 3;
    xv[1] = (xv[1] + xv[1] * xv[1]) / 2;

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // cell-baced velocities and Jacobian matrices
  // test piecewise linear deformation (part II)
  VectorPolynomial det;
  VectorPolynomial uc;
  MatrixPolynomial J;

  VectorPolynomial moments(2, 2);
  auto numi = std::make_shared<NumericalIntegration>(mesh0);
  std::vector<const char*> list = {"SerendipityPk"};
  
  for (auto name : list) {
    double fac(0.5), volume = mesh1->cell_volume(cell);
    for (int k = 1; k < 4; ++k) {
      // collect geometric data
      Teuchos::ParameterList plist;
      plist.set<std::string>("method", "Lagrange serendipity")
           .set<int>("method order", k)
           .set<std::string>("projector", "L2");
      auto maps = std::make_shared<MeshMaps_VEM>(mesh0, mesh1, plist);

      std::vector<WhetStone::VectorPolynomial> vf(nfaces);
      for (int f = 0; f < nfaces; ++f) maps->VelocityFace(f, vf[f]);

      // run differet discretization methods
      if (std::strcmp(name, "SerendipityPk") == 0) {
        maps->VelocityCell(cell, vf, vf, uc);
      }
      maps->Jacobian(uc, J);
      J(0, 0)(0) += 1.0;
      J(1, 1)(0) += 1.0;
      maps->Determinant(J, det);

      double tmp = numi->IntegratePolynomialCell(cell, det[0]);
      double err = tmp - volume;
      fac /= (k + 1);
      CHECK(fabs(err) < fac);

      uc.ChangeOrigin(mesh0->cell_centroid(cell));
      printf("k=%d  %14s  vol=%8.6g  err=%12.8f  |poly|=%9.6g %9.6g\n",
          k, name, tmp, err, uc[0].NormInf(), uc[1].NormInf());
    }
  }
}


/* ****************************************************************
* Comparison of least-square and Serendipity reconstructions
**************************************************************** */
TEST(DG_MAP_LEAST_SQUARE_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Comparison of reconstruction." << std::endl;
  auto comm = Amanzi::getDefaultComm();

  // create two meshes
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<Mesh> mesh0 = meshfactory.create("test/one_pentagon.exo");
  Teuchos::RCP<Mesh> mesh1 = meshfactory.create("test/one_pentagon.exo");

  // deform the second mesh
  int d(2), cell(0), nnodes(5), nfaces(5);
  AmanziGeometry::Point xv(d), yv(d);
  Entity_ID_List nodeids, faces;
  AmanziGeometry::Point_List new_positions, final_positions;

  // -- deformation function
  double dt(0.05);
  VectorPolynomial u(d, d); 

  for (int i = 0; i < d; ++i) {
    u[i].Reshape(d, 2, true);
    u[i](0, 0) = 1.0 + i;
    u[i](1, 0) = 2.0;
    u[i](1, 1) = 3.0 + i;
    u[i](2, 0) = 4.0;
    u[i](2, 1) = 5.0 + i;
    u[i](2, 2) = 6.0;
    // u[i](3, 1) = 6.0;
    // u[i](3, 3) = 7.0;
  }

  u *= dt;

  // -- deformation velocity on mesh faces
  std::vector<VectorPolynomial> vf(nfaces);
  VectorPolynomial vc1(d, d), vc2(d, d);

  for (int n = 0; n < nfaces; ++n) {
    vf[n] = u;
  }

  // -- mesh deformation
  for (int v = 0; v < nnodes; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    yv[0] = xv[0] + u[0].Value(xv);
    yv[1] = xv[1] + u[1].Value(xv);

    nodeids.push_back(v);
    new_positions.push_back(yv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // least-square calculation
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", "Lagrange serendipity")
       .set<int>("method order", 2)
       .set<std::string>("projector", "least square");
  auto maps = std::make_shared<MeshMaps_VEM>(mesh0, mesh1, plist);

  maps->VelocityCell(0, vf, vf, vc1);
  vc1.ChangeOrigin(mesh0->cell_centroid(cell));

  // Serendipity calculation
  plist.set<std::string>("method", "Lagrange serendipity")
       .set<int>("method order", 2)
       .set<std::string>("projector", "L2");
  maps = std::make_shared<MeshMaps_VEM>(mesh0, mesh1, plist);

  maps->VelocityCell(0, vf, vf, vc2);

  vc1 -= vc2;
  CHECK_CLOSE(0.0, vc1.NormInf(), 1e-12);
}


/* ****************************************************************
* Geometric conservation law (GCL) equation: dj/dt = j div u.
**************************************************************** */
TEST(DG_MAP_GCL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: GCL" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  // create 2D mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  // create a map
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", "Lagrange serendipity")
       .set<int>("method order", 2)
       .set<std::string>("projector", "L2");

  auto maps = std::make_shared<MeshMaps_VEM>(mesh, plist);

  // create a velocity polynomial for testing
  double t(0.5), dt(1e-2);
  VectorPolynomial u(2, 2, 3);
  for (int i = 0; i < 10; ++i) {
    u[0](i) = i;
    u[1](i) = (i * i) / 3 + 1.0;
  }

  // left-hand side of GCL, 2nd-order FD is exact
  VectorPolynomial det0, det1, det2;
  MatrixPolynomial J, Jt, C;

  maps->Jacobian(u, J);

  Jt = J * (t - 2 * dt);
  for (int i = 0; i < 2; ++i) Jt(i, i)(0) += 1.0; 
  maps->Determinant(Jt, det2);

  Jt = J * (t - dt);
  for (int i = 0; i < 2; ++i) Jt(i, i)(0) += 1.0; 
  maps->Determinant(Jt, det1);

  Jt = J * t;
  for (int i = 0; i < 2; ++i) Jt(i, i)(0) += 1.0; 
  maps->Determinant(Jt, det0);

  det0 = (1.5 * det0 - 2 * det1 + 0.5 * det2) * (1.0 / dt);

  // rigth-hand side of GCL
  VectorPolynomial v(u);

  maps->Cofactors(Jt, C);
  C.Multiply(u, v, true);
  det1 = Divergence(v);

  // error
  det1 -= det0;
  double err = det1.NormInf();
  std::cout << "GCL error=" << err << ", poly order=3" << std::endl;
  CHECK(err < 1e-10);

  // Piola compatibility condition
  for (int i = 0; i < 2; ++i) {
    v[0] = C(i, 0);
    v[1] = C(i, 1);
    u = Divergence(v);
    err = u.NormInf();
    std::cout << "Piola compatibility condition error=" << err << std::endl;
  }
}


/* ****************************************************************
* Hierachical velocity reconstruction in 3D
**************************************************************** */
TEST(DG_MAP_VELOCITY_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Velocity reconstruction in 3D." << std::endl;
  auto comm = Amanzi::getDefaultComm();

  // create two meshes
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<Mesh> mesh0 = meshfactory.create("test/cube_unit.exo", true, true);
  Teuchos::RCP<Mesh> mesh1 = meshfactory.create("test/cube_unit.exo", true, true);

  // deform the second mesh
  int d(3), nnodes(8), nfaces(6), nedges(12);
  AmanziGeometry::Point xv(d), yv(d);
  Entity_ID_List nodeids, edges, faces;
  AmanziGeometry::Point_List new_positions, final_positions;

  // -- deformation function
  int order(1);
  VectorPolynomial u(d, d); 
  for (int i = 0; i < d; ++i) {
    u[i].Reshape(d, order, true);
    u[i](0, 0) = 1.0 + i;
    u[i](1, 0) = 2.0 - i;
    u[i](1, 1) = 3.0;
    u[i](1, 2) = 4.0 + 2 * i;
    u[i].set_origin(AmanziGeometry::Point(d));
  }
  u *= 0.05;

  for (int v = 0; v < nnodes; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    for (int i = 0; i < d; ++i) 
      yv[i] = xv[i] + u[i].Value(xv);

    nodeids.push_back(v);
    new_positions.push_back(yv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // velocity calculation
  // -- on edges
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", "Lagrange serendipity")
       .set<int>("method order", order)
       .set<std::string>("projector", "H1");
  auto maps = std::make_shared<MeshMaps_VEM>(mesh0, mesh1, plist);

  std::vector<VectorPolynomial> ve(nedges); 
  mesh0->cell_get_edges(0, &edges);
  for (int n = 0; n < nedges; ++n) {
    maps->VelocityEdge(edges[n], ve[n]);
  }

  // -- on faces
  std::vector<VectorPolynomial> vf(nfaces); 
  mesh0->cell_get_faces(0, &faces);
  for (int n = 0; n < nfaces; ++n) {
    maps->VelocityFace(faces[n], vf[n]);
  }

  // -- in cell
  VectorPolynomial vc;
  maps->VelocityCell(0, ve, vf, vc);
  vc.ChangeOrigin(AmanziGeometry::Point(3));

  vc -= u;
  CHECK_CLOSE(0.0, vc.NormInf(), 1e-12);
}


