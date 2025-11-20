/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Operators

*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "KDTree.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshVirtual.hh"
#include "Op_Face_Cell.hh"
#include "Operator_Cell.hh"
#include "PDE_Diffusion.hh"
#include "PDE_ElasticityCurvedFace.hh"
#include "Tensor.hh"

// Operators
#include "AnalyticElasticity01.hh"
#include "AnalyticElasticity03.hh"

#include "operator_virtual.hh"

using namespace Amanzi;
using namespace Teuchos;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;


/* *****************************************************************
* Test
***************************************************************** */
double
RunTestVirtualElasticity(int d, int level = 1)
{
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: elasticity solver on a virtual logical mesh \n";

  // read parameter list
  std::string xmlFileName = "test/operator_virtual_elasticity.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // read input parameters
  int nx, ny, nz;
  double Lx, Ly, Lz;
  nx = level * plist->get<int>("mesh cells in x direction", 20);
  ny = level * plist->get<int>("mesh cells in y direction", 20);
  Lx = plist->get<double>("domain size in x direction", 1.0);
  Ly = plist->get<double>("domain size in y direction", 1.0);
  std::vector<double> origin = plist->get<Teuchos::Array<double>>("origin").toVector();

  if (d == 3) {
    nz = level * plist->get<int>("mesh cells in z direction", 20);
    Lz = plist->get<double>("domain size in y direction", 1.0);
    origin.push_back(1.0);
  }

  // generate cloud of points
  std::vector<double> cell_volumes;
  Point_List cell_centroids;
  Point_List face_centroids_bnd;
  Point_List face_normals_bnd;

  if (d == 2) {
    generatePointCloud2D(nx, ny,
                         Lx, Ly, 
                         origin,
                         face_centroids_bnd, 
                         face_normals_bnd, 
                         cell_centroids,
                         cell_volumes);
  } else {
    generatePointCloud3D(nx, ny, nz,
                         Lx, Ly, Lz,
                         origin,
                         face_centroids_bnd, 
                         face_normals_bnd, 
                         cell_centroids,
                         cell_volumes);
  }

  auto face_cells = createConnectivityGraph(face_centroids_bnd, cell_centroids);

  int ncells = cell_centroids.size();
  int nfaces = face_cells.size();
  int nfaces_bnd = face_centroids_bnd.size();
  int nfaces_int = nfaces - nfaces_bnd;

  // create first virtual mesh for geometry moments
  std::vector<AmanziGeometry::Point> face_centroids(nfaces), face_normals(nfaces);
  for (int f = 0; f < nfaces; ++f) {
    int c1 = face_cells[f][0];
    if (face_cells[f].size() == 1) {
      int fb = f - nfaces_int;
      face_centroids[f] = face_centroids_bnd[fb];
      face_normals[f] = face_normals_bnd[fb];
    } else {
      int c2 = face_cells[f][1];
      face_centroids[f] = (cell_centroids[c1] + cell_centroids[c2]) / 2;
      face_normals[f] = cell_centroids[c2] - cell_centroids[c1];
    }
  }

  auto mesh_fw1 = Teuchos::rcp(new MeshVirtual(comm,
                                               plist,
                                               face_cells,
                                               cell_centroids,
                                               cell_volumes,
                                               face_centroids,
                                               face_normals));

  auto mesh1 = Teuchos::rcp(new AmanziMesh::MeshCache<MemSpace_kind::HOST>(
    mesh_fw1, Teuchos::rcp(new MeshVirtualAlgorithms()), Teuchos::null));

  auto moments = computeGeometryMoments(mesh1, nfaces_int, nfaces, plist);
  auto& moments_f = *moments->ViewComponent("face");

  // -- verification
  AmanziGeometry::Point normal(d);
  for (int c = 0; c < ncells; ++c) {
    const auto& [faces, dirs] = mesh1->getCellFacesAndDirections(c);
    int nfaces = faces.size();

    double sum(0.0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      for (int k = 0; k < d; ++k) normal[k] = moments_f[k][f];
      sum += normal[0] * dirs[n];
    }
    CHECK_CLOSE(0.0, sum, 1e-8);
  }

  // create second virtual mesh with fully complete geometry
  for (int f = 0; f < nfaces; ++f) {
    for (int k = 0; k < d; ++k) face_normals[f][k] = moments_f[k][f];
  }

  auto mesh_fw2 = Teuchos::rcp(new MeshVirtual(comm,
                                               plist,
                                               face_cells,
                                               cell_centroids,
                                               cell_volumes,
                                               face_centroids,
                                               face_normals));

  auto mesh2 = Teuchos::rcp(new AmanziMesh::MeshCache<MemSpace_kind::HOST>(
    mesh_fw2, Teuchos::rcp(new MeshVirtualAlgorithms()), Teuchos::null));

  // -- initialize diffusion problem
  Teuchos::ParameterList pde_list;
  pde_list.set<double>("diagonal shift", 0.0);
  pde_list.set<std::string>("matrix type", "stiffness");
  pde_list.sublist("schema")
    .set<std::string>("method", "elasticity weak symmetry curved face")
    .set<int>("method order", 1)
    .set<std::string>("base", "cell");
  auto pde = Teuchos::rcp(new Operators::PDE_ElasticityCurvedFace(pde_list, mesh2));
  auto op = pde->global_operator();

  // -- set diffusion coefficient
  auto ana = Teuchos::RCP(new AnalyticElasticity01(mesh2, 1.0, 1.0));
  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < cell_centroids.size(); ++c) {
    const AmanziGeometry::Point& xc = mesh2->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana->Tensor(xc, 0.0);
    K->push_back(Kc);
  }
  pde->Setup(K, false);

  // -- verification
  auto bf = pde->get_bf();
  for (int c = 0; c < ncells; ++c) {
    const auto& [faces, dirs] = mesh2->getCellFacesAndDirections(c);
    int nfaces = faces.size();

    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        double sum(0.0);
        for (int n = 0; n < nfaces; ++n) {
          int f = faces[n];
          const auto& normal = mesh2->getFaceNormal(f);
          const AmanziGeometry::Point& xf = (*bf)[f];
          sum += normal[i] * dirs[n] * xf[j];
        }
        CHECK_CLOSE(sum, ((i == j) ? mesh2->getCellVolume(c) : 0.0), 1e-10);
      }
    }
  }

  // -- finalize diffusion problem
  auto bc = Teuchos::rcp(new Operators::BCs(mesh2, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
  auto& bc_model = bc->bc_model();
  auto& bc_value = bc->bc_value_vector(d);

  for (int f = nfaces_int; f < nfaces; ++f) {
    const auto& xf = (*bf)[f];
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;

    auto tmp = ana->velocity_exact(xf, 0.0);
    for (int k = 0; k < d; ++k) bc_value[f][k] = tmp[k];
  }
  pde->AddBCs(bc, bc);

  op->Init();

  auto rhs = op->rhs();
  auto& rhs_c = *rhs->ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
    const auto& xc = mesh2->getCellCentroid(c);
    const auto& vol = mesh2->getCellVolume(c);
    const auto& tmp = ana->source_exact(xc, 0.0);
    for (int k = 0; k < d; ++k) rhs_c[k][c] = vol * tmp[k];
  }

  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->ApplyBCs(true, true, true);

  op->SymbolicAssembleMatrix();
  op->AssembleMatrix();
  // std::cout << *op->A() << std::endl; exit(0);

  // -- solver
  auto pc_list = Teuchos::sublist(plist, "preconditioners", true);
  auto solver_list = Teuchos::sublist(plist, "solvers", true);
  op->set_inverse_parameters("PCG", *solver_list, "Hypre AMG", *pc_list);

  auto sol(*rhs);
  op->ApplyInverse(*rhs, sol);

  double l2_err(0.0), inf_err(0.0), unorm(0.0);
  auto& sol_c = *sol.ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    const auto& xc = mesh2->getCellCentroid(c);
    double vol = mesh2->getCellVolume(c);
    auto tmp = ana->velocity_exact(xc, 0.0);

    for (int k = 0; k < d; ++k) {
      inf_err = std::max(inf_err, std::fabs(tmp[k] - sol_c[k][c]));
      l2_err += vol * std::pow(tmp[k] - sol_c[k][c], 2);
      unorm += vol * std::pow(sol_c[k][c], 2);
    }
  }
  printf("L2(p) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", std::sqrt(l2_err / unorm), inf_err, std::sqrt(unorm));
  AMANZI_ASSERT(l2_err < 1e-8);

  // -- visualization
  std::ofstream file1("xyc");
  for (int c = 0; c < ncells; ++c) {
    file1 << mesh2->getCellCentroid(c);
    for (int k = 0; k < d; ++k) file1 << " " << sol_c[k][c];
    file1 << std::endl;
  }
  file1.close();

  std::ofstream file2("xyn");
  for (int f = 0; f < nfaces; ++f) {
    file2 << (*bf)[f] << " " << mesh2->getFaceNormal(f) << std::endl;
  }
  file2.close();

  return l2_err;
}


TEST(OPERATOR_VIRTUAL_ELASTICITY_2D)
{
  RunTestVirtualElasticity(2, 1);
}


TEST(OPERATOR_VIRTUAL_ELASTICITY_3D)
{
  RunTestVirtualElasticity(3);
}
