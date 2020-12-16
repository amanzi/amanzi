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
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "SingleFaceMesh.hh"
#include "SurfaceCoordinateSystem.hh"
#include "VectorObjects.hh"
#include "VEM_NedelecSerendipity.hh"
#include "VEM_RaviartThomasSerendipity.hh"

// Amanzi::Operators
#include "AnalyticElectromagnetics05.hh"
#include "AnalyticElectromagnetics02.hh"
#include "MeshDeformation.hh"


/* *****************************************************************
* Test approximation of || P1_c(E) - E ||.
* **************************************************************** */
void ProjectorRTAccuracy() {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: global accuracy of the L2 RT projector" << std::endl;

  // create a MSTK mesh framework
  ParameterList plist;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, plist, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  int order(1);
  int ndf = PolynomialSpaceDimension(2, order);

  WhetStone::Tensor K(3, 1);
  K(0, 0) = 1.0;

  for (int nx = 1; nx < 3; nx *= 2) {
    // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, true, true);
    // DeformMesh(mesh, 5, 0.1);
    std::string name = (nx == 1) ? "4" : "8";
    RCP<Mesh> mesh = meshfactory.create("test/hexes" + name + ".exo", true, true);
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    AnalyticElectromagnetics02 ana(1.0, mesh);
    WhetStone::NumericalIntegration numi(mesh);

    plist.set<int>("method order", order);
    WhetStone::VEM_RaviartThomasSerendipity vem(plist, mesh);

    double err1(0.0), err2(0.0), err3(0.0), err4(0.0);
    WhetStone::DenseMatrix A;
    std::vector<int> dirs;
    AmanziMesh::Entity_ID_List faces, nodes;
    AmanziGeometry::Point zero(2);

    for (int c = 0; c < ncells_owned; ++c) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);

      // compute DOFs
      mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      std::vector<double> moments;
      DenseVector v1(nfaces * ndf), v2(nfaces * ndf), v3(nfaces * ndf);
      std::vector<VectorPolynomial> vf(nfaces);

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        double area = mesh->face_area(f);
        const auto& xf = mesh->face_centroid(f);
        const auto& normal = mesh->face_normal(f);
        auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);

        ana.set_parameters(normal / area, 2, 0.0);
        numi.CalculateFunctionMomentsFace(f, &ana, order, moments, 4);
        for (int k = 0; k < moments.size(); ++k) {
          v1(ndf * n + k) = moments[k] * dirs[n];
        }

        vf[n].Reshape(2, 3, order, true);
        vf[n].set_origin(zero);

        std::vector<double> scale;
        Polynomial poly(2, 0);
        poly(0) = 1.0;
        numi.CalculatePolynomialMomentsFace(f, poly, 2 * order, scale);

        DenseMatrix M(3, 3);
        M.PutScalar(0.0);

        for (int i = 0; i < 3; ++ i) {
          M(0, 0) = scale[0];
          M(1, 1) = scale[3];
          M(1, 2) = M(2, 1) = scale[4];
          M(2, 2) = scale[5];
          M.Inverse();

          double f1 = normal[i] / area;
          double f2 = f1 / std::sqrt(area);
          vf[n][i](0) = f1 * M(0, 0) * moments[0];
          vf[n][i](1) = f2 * (M(1, 1) * moments[1] + M(1, 2) * moments[2]);
          vf[n][i](2) = f2 * (M(2, 1) * moments[1] + M(2, 2) * moments[2]);

          vf[n][i].InverseChangeCoordinates(xf, *coordsys->tau());  
        }
      }

      // compute projector
      VectorPolynomial Bc;
      vem.L2Cell(c, vf, vf, nullptr, Bc);

      // error || Bc(xc) - B(xc) ||
      auto B1 = ana.magnetic_exact(xc, 0.0);
      auto B2 = Bc.Value(xc);
      for (int i = 0; i < 3; ++i) {
        err1 += std::pow(B1[i] - B2(i), 2.0) * mesh->cell_volume(c);
      }


      // compute DOFs of the projector for numerical integration
      for (int n = 0; n < nfaces; ++n) vf[n] = Bc;
      vem.CalculateDOFsOnBoundary(c, vf, vf, v2);
      v3 = v2;

      // error ||| Ec - E |||_E
      DenseMatrix Mf;
      vem.MassMatrix(c, K, Mf);
      v1 -= v2;
      Mf.Multiply(v1, v2, false);
      err2 += v1 * v2;

      // error || Ec - E ||
      for (int i = 0; i < 3; ++i) {
        err3 += std::pow(v1(i), 2.0) * mesh->cell_volume(c);
      }


      // orthogonality v1 . P1(Bc)^I
      err4 += fabs(v1 * v3) * mesh->cell_volume(c);
    }
    std::cout << "nx=" << nx << " errors=" << std::sqrt(err1) << " " 
              << std::sqrt(err2) << " " << std::sqrt(err3) 
              << "  ortho=" << err4 << std::endl;
    CHECK_CLOSE(err4, 0.0, 1e-10);
  }
}


TEST(PROJECTOR_RT_ACCURACY) {
  ProjectorRTAccuracy();
}


/* *****************************************************************
* Test approximation of || P1_c(E) - E ||.
* **************************************************************** */
void ProjectorNDAccuracy() {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: global accuracy of the ND L2 projector" << std::endl;

  // create a MSTK mesh framework
  ParameterList plist;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, plist, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  int order(1);
  int nde = PolynomialSpaceDimension(1, order);

  WhetStone::Tensor K(3, 1);
  K(0, 0) = 1.0;

  for (int nx = 1; nx < 3; nx *= 2) {
    // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, true, true);
    // DeformMesh(mesh, 5, 0.1);
    std::string name = (nx == 1) ? "4" : "8";
    RCP<Mesh> mesh = meshfactory.create("test/hexes" + name + ".exo", true, true);
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    AnalyticElectromagnetics02 ana(1.0, mesh);
    WhetStone::NumericalIntegration numi(mesh);

    plist.set<int>("method order", order)
         .set<int>("type", 2);
    WhetStone::VEM_NedelecSerendipity vem(plist, mesh);

    double err1(0.0), err2(0.0), err4(0.0);
    WhetStone::DenseMatrix A;
    AmanziMesh::Entity_ID_List edges;
    AmanziGeometry::Point zero(1);

    for (int c = 0; c < ncells_owned; ++c) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);

      // compute DOFs
      mesh->cell_get_edges(c, &edges);
      int nedges = edges.size();

      std::vector<double> moments;
      DenseVector v1(nedges * nde), v2(nedges * nde), v3(nedges * nde);
      std::vector<VectorPolynomial> ve(nedges);

      for (int n = 0; n < nedges; ++n) {
        int e = edges[n];
        double len = mesh->edge_length(e);
        const auto& xe = mesh->edge_centroid(e);
        const auto& tau = mesh->edge_vector(e);
        std::vector<AmanziGeometry::Point> vtau(1, tau);

        ana.set_parameters(tau / len, 0, 0.0);
        numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 4);
        for (int k = 0; k < moments.size(); ++k) {
          v1(nde * n + k) = moments[k];
        }

        ve[n].Reshape(1, 3, order, true);
        ve[n].set_origin(zero);

        for (int i = 0; i < 3; ++ i) {
          for (int k = 0; k < moments.size(); ++k) {
            double l22 = (k == 0) ? 1.0 : 1.0 / 12;
            ve[n][i](k) = moments[k] * tau[i] / len / l22;
          }
          ve[n][i].InverseChangeCoordinates(xe, vtau);  
        }
      }

      // compute projector
      VectorPolynomial Ec;
      vem.L2Cell(c, ve, ve, nullptr, Ec);

      // error || Ec(xc) - E(xc) ||
      auto E1 = ana.electric_exact(xc, 0.0);
      auto E2 = Ec.Value(xc);
      for (int i = 0; i < 3; ++i) {
        err1 += std::pow(E1[i] - E2(i), 2.0) * mesh->cell_volume(c);
      }


      // compute DOFs of the projector for numerical integration
      for (int n = 0; n < nedges; ++n) ve[n] = Ec;
      vem.CalculateDOFsOnBoundary(mesh, c, ve, ve, v2);
      v3 = v2;

      // error ||| Ec - E |||_E
      DenseMatrix Me;
      vem.MassMatrix(c, K, Me);
      v1 -= v2;
      Me.Multiply(v1, v2, false);
      err2 += v1 * v2;


      // orthogonality v1 . P1(Ec)^I
      err4 += fabs(v1 * v3) * mesh->cell_volume(c);
    }
    std::cout << "nx=" << nx << " errors=" << std::sqrt(err1) << " " 
              << std::sqrt(err2) << "  ortho=" << err4 << std::endl;

    CHECK_CLOSE(err4, 0.0, 1e-10);
  }
}


TEST(PROJECTOR_ND_ACCURACY) {
  ProjectorNDAccuracy();
}


/* *****************************************************************
* Test approximation of B = curl E with the primary curl operator 
* CURL is exact, so error depends only on accuracy of a quadrature.
* **************************************************************** */
void PrimaryCurl() {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: PRIMARY operator curl" << std::endl;

  // create a MSTK mesh framework
  ParameterList plist;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, plist, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  int order(1), type(1);
  for (int nx = 4; nx < 5; nx *= 2) {
    // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, true, true);
    // DeformMesh(mesh, 1, 0.1);
    RCP<Mesh> mesh = meshfactory.create("test/hexes4.exo", true, true);
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    AnalyticElectromagnetics02 ana(1.0, mesh);
    WhetStone::NumericalIntegration numi(mesh);

    plist.set<int>("method order", order)
         .set<int>("type", type);
    WhetStone::VEM_NedelecSerendipity vem(plist, mesh);

    double err(0.0);
    WhetStone::DenseMatrix A;
    std::vector<int> fdirs;
    AmanziMesh::Entity_ID_List faces, edges;

    for (int c = 0; c < ncells_owned; ++c) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);

      vem.CurlMatrix(c, A);

      // electric field moments (on edges)
      std::vector<double> moments, Eex, Bex;
      mesh->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n < nedges; ++n) {
        int e = edges[n];
        double len = mesh->edge_length(e);
        const auto& tau = mesh->edge_vector(e);

        ana.set_parameters(tau / len, 0, 0.0);
        numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 4);
        for (int k = 0; k < moments.size(); ++k) Eex.push_back(moments[k]);
      }

      // electric field moments (on faces)
      mesh->cell_get_faces_and_dirs(c, &faces, &fdirs);
      int nfaces = faces.size();

      if (type == 1) {
        WhetStone::Polynomial pf(2, order);

        for (int n = 0; n < nfaces; ++n) {
          int f = faces[n];
          double area = mesh->face_area(f);
          const auto& xf = mesh->face_centroid(f);
          const auto& normal = mesh->face_normal(f);
          WhetStone::SurfaceCoordinateSystem coordsys(xf, normal);

          for (auto it = pf.begin(1); it < pf.end(); ++it) {
            double factor = std::pow(area, -(double)it.MonomialSetOrder() / 2);
            WhetStone::Polynomial fmono(2, it.multi_index(), factor);
            auto rot = Rot2D(fmono);

            AmanziGeometry::Point tmp = rot[0](0) * (*coordsys.tau())[0] 
                                      + rot[1](0) * (*coordsys.tau())[1];
            ana.set_parameters(tmp, 0, 0.0);
            numi.CalculateFunctionMomentsFace(f, &ana, order - 1, moments, 4);
            for (int k = 0; k < moments.size(); ++k) Eex.push_back(moments[k]);
          }
        }
      }
      
      // magnetic field moments wrt exterior normal
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        double area = mesh->face_area(f);
        AmanziGeometry::Point normal = mesh->face_normal(f);

        ana.set_parameters(normal / area, 2, 0.0);
        numi.CalculateFunctionMomentsFace(f, &ana, order, moments, 4);
        for (int k = 0; k < moments.size(); ++k) {
          Bex.push_back(moments[k] * fdirs[n]);
        }
      }

      // error
      WhetStone::DenseVector Eh(Eex), Bh(Bex);
      A.Multiply(Eh, Bh, false);
      for (int i = 0; i < Bh.NumRows(); ++i) {
        err += std::pow(Bh(i) - Bex[i], 2.0) * mesh->cell_volume(c);
      }
    }
    std::cout << "nx=" << nx << " " << std::sqrt(err) << std::endl;
    CHECK_CLOSE(err, 0.0, 1e-10);
  }
}


TEST(PRIMARY_OPERATOR_CURL) {
  PrimaryCurl();
}


/* *****************************************************************
* Test approximation of E = curl B with the dual curl operator 
* **************************************************************** */
void DualCurl() {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::WhetStone;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: DUAL operator curl" << std::endl;

  // create a MSTK mesh framework
  ParameterList plist;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, plist, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  int order(1);
  int nde = PolynomialSpaceDimension(1, order);

  for (int nx = 1; nx < 2; nx *= 2) {
    RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, true, true);
    // DeformMesh(mesh, 1, 0.1);
    // RCP<Mesh> mesh = meshfactory.create("test/one_trapezoid.exo", true, true);
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    AnalyticElectromagnetics05 ana(mesh);
    WhetStone::NumericalIntegration numi(mesh);

    plist.set<int>("method order", order)
         .set<int>("type", 2);
    WhetStone::VEM_NedelecSerendipity vem_nd(plist, mesh);
    WhetStone::VEM_RaviartThomasSerendipity vem_rt(plist, mesh);

    double err(0.0);
    WhetStone::DenseMatrix C, Me, Mf;
    std::vector<int> fdirs, edirs, map;
    AmanziMesh::Entity_ID_List faces, edges, fedges;

    WhetStone::Tensor K(3, 1);
    K(0, 0) = 1.0;

    for (int c = 0; c < ncells_owned; ++c) {
      vem_nd.CurlMatrix(c, C);
      vem_nd.MassMatrix(c, K, Me);
      vem_rt.MassMatrix(c, K, Mf);

      // electric field moments
      std::vector<double> moments, Eex, Bex;
      mesh->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n < nedges; ++n) {
        int e = edges[n];
        double len = mesh->edge_length(e);
        const auto& tau = mesh->edge_vector(e);

        ana.set_parameters(tau / len, 0, 0.0);
        numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 4);
        for (int k = 0; k < moments.size(); ++k) Eex.push_back(moments[k]);
      }

      // magnetic field moments wrt exterior normal
      mesh->cell_get_faces_and_dirs(c, &faces, &fdirs);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        double area = mesh->face_area(f);
        AmanziGeometry::Point normal = mesh->face_normal(f);

        ana.set_parameters(normal / area, 2, 0.0);
        numi.CalculateFunctionMomentsFace(f, &ana, order, moments, 4);
        for (int k = 0; k < moments.size(); ++k) Bex.push_back(moments[k] * fdirs[n]);
      }

      // Since the dual mimetic operator uses inverse mass matrix, we test
      // approximation using the discrete equation Me E = CURL^T Mf B.
      WhetStone::DenseVector E1(Eex), E2(Eex), E3(Eex), B1(Bex), B2(Bex);
      Me.Multiply(E1, E2, false);

      Mf.Multiply(B1, B2, false);
      C.Multiply(B2, E3, true);

      // Boundary terms
      DenseMatrix Se;
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& xf = mesh->face_centroid(f);
        AmanziGeometry::Point normal = mesh->face_normal(f);
        double area = mesh->face_area(f);

        if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
            fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6 ||
            fabs(xf[2]) < 1e-6 || fabs(1.0 - xf[2]) < 1e-6 || true) {

          vem_nd.MassMatrixFace(f, K, Se);

          mesh->face_to_cell_edge_map(f, c, &map);
          mesh->face_get_edges_and_dirs(f, &fedges, &edirs);
          int nfedges = fedges.size();

          DenseVector v1(nfedges * nde), v2(nfedges * nde);

          for (int i = 0; i < nfedges; ++i) {
            int e = fedges[i];
            const AmanziGeometry::Point& tau = mesh->edge_vector(e);
            double len = mesh->edge_length(e);

            ana.set_parameters((normal ^ tau) / area / len, 2, 0.0);
            numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 4);
            for (int k = 0; k < nde; ++k) {
              v1(i * nde + k) = moments[k] * fdirs[n];
            }
          }

          Se.Multiply(v1, v2, false);

          for (int i = 0; i < nfedges; ++i) {
            for (int k = 0; k < nde; ++k) {
              E3(map[i] * nde + k) -= v2(i * nde + k); 
            }
          }
        }
      }

      for (int i = 0; i < E3.NumRows(); ++i) {
        err += std::pow(E2(i) - E3(i), 2.0) / mesh->cell_volume(c);
      }
    }
    std::cout << "nx=" << nx << " " << std::sqrt(err) << std::endl;
  }
}


TEST(DUAL_OPERATOR_CURL) {
  DualCurl();
}


