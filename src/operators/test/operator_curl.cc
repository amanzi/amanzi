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
#include "VEM_NedelecSerendipityType2.hh"

// Amanzi::Operators
#include "AnalyticElectromagnetics01.hh"
#include "AnalyticElectromagnetics02.hh"
#include "MeshDeformation.hh"

/* *****************************************************************
* TBW 
* **************************************************************** */
void Curl2ndOrder() {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2nd order operator curl" << std::endl;

  // create a MSTK mesh framework
  ParameterList plist;
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, plist, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  int order(1);
  for (int nx = 2; nx < 3; nx *= 2) {
    // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, true, true);
    // DeformMesh(mesh, 1, 0.1);
    RCP<Mesh> mesh = meshfactory.create("test/kershaw08.exo", true, true);
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    AnalyticElectromagnetics02 ana(1.0, mesh);
    WhetStone::NumericalIntegration<AmanziMesh::Mesh> numi(mesh);

    plist.set<int>("method order", order);
    WhetStone::VEM_NedelecSerendipityType2 vem(plist, mesh);

    double err(0.0);
    WhetStone::DenseMatrix A;
    std::vector<int> fdirs;
    AmanziMesh::Entity_ID_List faces, edges;

    for (int c = 0; c < ncells_owned; ++c) {
      vem.CurlMatrix(c, A);

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

if(f==49){
std::vector<int> edirs;
mesh->face_get_edges_and_dirs(f, &edges, &edirs);
for(int i = 0; i < edges.size(); ++i) {
  int e = edges[i];
  double len = mesh->edge_length(e);
  const auto& tau = mesh->edge_vector(e);

  ana.set_parameters(tau / len, 0, 0.0);
  numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 4);
std::cout << moments[0] << "\n " << moments[1] << std::endl;
}
}

      }

      // error
      WhetStone::DenseVector Eh(Eex), Bh(Bex);
      A.Multiply(Eh, Bh, false);
      // for (int i = 0; i < Bh.NumRows(); ++i) {
      for (int i = 2; i < Bh.NumRows(); i += 3) {
        err += std::pow(Bh(i) - Bex[i], 2.0) * mesh->cell_volume(c);
// std::cout << c << " " << i << " " << Bh(i) << " " << Bex[i] << " " << Bh(i) / (fabs(Bex[i]) + 1e-15) << std::endl;
      }
    }
    std::cout << "nx=" << nx << " " << std::sqrt(err) << std::endl;
  }
}


TEST(OPERATOR_CURL_2ND_ORDER) {
  Curl2ndOrder();
}

