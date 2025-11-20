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
#include "PDE_AdvectionUpwind.hh"
#include "PDE_DiffusionCurvedFace.hh"
#include "Tensor.hh"

// Operators
#include "Analytic00.hh"
#include "Analytic03.hh"

#include "operator_virtual.hh"

using namespace Amanzi;
using namespace Teuchos;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;


/* *****************************************************************
* Test
***************************************************************** */
void
RunTestVirtualAdvectionDiffusion(int d)
{
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: elliptic solver on a virtual logical mesh \n";

  // read parameter list
  std::string xmlFileName = "test/operator_virtual_advdiff.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // read input parameters
  int nx, ny, nz;
  double Lx, Ly, Lz;
  nx = plist->get<int>("mesh cells in x direction", 20);
  ny = plist->get<int>("mesh cells in y direction", 20);
  Lx = plist->get<double>("domain size in x direction", 1.0);
  Ly = plist->get<double>("domain size in y direction", 1.0);
  std::vector<double> origin = plist->get<Teuchos::Array<double>>("origin").toVector();

  if (d == 3) {
    nz = plist->get<int>("mesh cells in z direction", 20);
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
  AmanziGeometry::Point vel(1.0, 2.0);
  auto ana = Teuchos::RCP(new Analytic00(mesh2, 2, 0.0, vel));

  Teuchos::ParameterList pde_list;
  pde_list.set<double>("penalty", 0.0);
  auto pde = Teuchos::rcp(new Operators::PDE_DiffusionCurvedFace(pde_list, mesh2, nullptr));
  auto op = pde->global_operator();

  // -- set diffusion coefficient
  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < cell_centroids.size(); ++c) {
    const AmanziGeometry::Point& xc = mesh2->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana->TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }
  pde->SetTensorCoefficient(K);

  // -- finalize diffusion problem
  auto bc = Teuchos::rcp(new Operators::BCs(mesh2, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  auto& bc_value = bc->bc_value();
  auto& bc_model = bc->bc_model();

  int dir;
  auto bf = pde->get_bf();
  for (int f = nfaces_int; f < nfaces; ++f) {
    const auto& xf0 = mesh2->getFaceCentroid(f); // original centroid
    const auto& xf = (*bf)[f];
    if (fabs(xf0[0] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana->pressure_exact(xf, 0.0);
    } else if (fabs(xf0[1] - 1.0) < 1e-6) { // inflow bottom boundary
      double area = mesh2->getFaceArea(f);
      auto normal = getFaceNormalExterior(*mesh2, f, &dir) / area;

      bc_model[f] = OPERATOR_BC_TOTAL_FLUX;
      double flux1 = ana->velocity_exact(xf, 0.0) * normal;
      double flux2 = (ana->advection_exact(xf, 0.0) * normal) * ana->pressure_exact(xf, 0.0);
      bc_value[f] = flux1 + flux2;
    } else if (fabs(xf0[0] - 2.0) < 1e-6 || fabs(xf0[1] - 2.0) < 1e-6) { // outflow right boundary
      double area = mesh2->getFaceArea(f);
      auto normal = getFaceNormalExterior(*mesh2, f, &dir);

      bc_model[f] = OPERATOR_BC_NEUMANN;
      bc_value[f] = (ana->velocity_exact(xf, 0.0) * normal) / area;
    } else {
      AMANZI_ASSERT(false);
    }
  }
  pde->AddBCs(bc, bc);

  auto rhs = op->rhs();
  auto& rhs_c = *rhs->ViewComponent("cell");

  for (int c = 0; c < ncells; ++c) {
    const auto& xc = mesh2->getCellCentroid(c);
    const auto& vol = mesh2->getCellVolume(c);
    rhs_c[0][c] = vol * ana->source_exact(xc, 0.0);
  }

  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->ApplyBCs(true, true, true);

  // -- finalize advection operator
  pde_list = plist->sublist("PK operator").sublist("advection operator");
  auto pde_adv = Teuchos::rcp(new PDE_AdvectionUpwind(pde_list, op));
  pde_adv->SetBCs(bc, bc);

  auto u = Teuchos::rcp(new CompositeVector(*rhs));
  auto& u_f = *u->ViewComponent("face");
  for (int f = 0; f < nfaces; f++) {
    auto tmp = ana->advection_exact((*bf)[f], 0.0);
    u_f[0][f] = tmp * mesh2->getFaceNormal(f);
  }

  pde_adv->Setup(*u);
  pde_adv->UpdateMatrices(u.ptr());
  pde_adv->ApplyBCs(false, true, false);

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
    double tmp = ana->pressure_exact(xc, 0.0);

    inf_err = std::max(inf_err, std::fabs(tmp - sol_c[0][c]));
    l2_err += vol * std::pow(tmp - sol_c[0][c], 2);
    unorm += vol * std::pow(tmp, 2);
  }
  printf("L2(p) =%10.7f  Inf =%9.6f  |u|=%9.6f\n", std::sqrt(l2_err / unorm), inf_err, std::sqrt(unorm));
  AMANZI_ASSERT(l2_err < 1e-3);

  // -- visualization
  std::ofstream file1("xyc");
  for (int c = 0; c < ncells; ++c) {
    file1 << mesh2->getCellCentroid(c) << " " << sol_c[0][c] << std::endl;
  }
  file1.close();
}


TEST(OPERATOR_VIRTUAL_ADVECTION_DIFFUSION)
{
  RunTestVirtualAdvectionDiffusion(2);
}

