/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "DG_Modal.hh"
#include "Explicit_TI_RK.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "OutputXDMF.hh"

// Amanzi::Operators
#include "RemapDG_Tests.hh"

#include "AnalyticDG04.hh"

namespace Amanzi {

class MyRemapDG : public RemapDG_Tests<AnalyticDG04> {
 public:
  MyRemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
            const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
            Teuchos::ParameterList& plist)
    : RemapDG_Tests<AnalyticDG04>(mesh0, mesh1, plist){};
  ~MyRemapDG(){};

  // access
  const std::vector<WhetStone::VectorPolynomial> det() const { return *det_; }
  const std::shared_ptr<WhetStone::MeshMaps> maps() const { return maps_; }
};

} // namespace Amanzi


/* *****************************************************************
 * Remap of polynomilas in two dimensions. Explicit time scheme.
 * Dual formulation places gradient and jumps on a test function.
 ***************************************************************** */
void
RemapTestsDualRK(const Amanzi::Explicit_TI::method_t& rk_method,
                 std::string map_name, std::string file_name, int nx, int ny,
                 int nz, double dt, int deform = 1)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  int dim = (nz == 0) ? 2 : 3;

  auto comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();

  // read parameter list
  std::string xmlFileName = "test/operator_remap.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  int order = plist.sublist("PK operator")
                .sublist("flux operator")
                .get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", map_name);

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  int vel_order = map_list.get<int>("method order");

  if (getRank == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string map_name = map_list.get<std::string>("map name");

    std::cout << "\nTest: " << dim << "D remap:"
              << " mesh=" << ((ny == 0) ? file_name : "square")
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order << ", map=" << map_name
              << std::endl;

    std::cout << "      map details: order=" << vel_order
              << ", projector=" << vel_projector << ", method=\"" << vel_method
              << "\"" << std::endl;
  }

  // create two meshes
  MeshFactory meshfactory(comm);
  // meshfactory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  meshfactory.set_preference(Preference({ AmanziMesh::Framework::MSTK }));

  Teuchos::RCP<const Mesh> mesh0;
  Teuchos::RCP<Mesh> mesh1;
  if (dim == 2 && ny != 0) {
    mesh0 = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, ny);
    mesh1 = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, ny);
  } else if (dim == 2) {
    mesh0 = meshfactory.create(file_name);
    mesh1 = meshfactory.create(file_name);
  } else {
    mesh0 =
      meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, true, true);
    mesh1 =
      meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, true, true);
  }

  int ncells_owned =
    mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // create and initialize cell-based field
  CompositeVectorSpace cvs1, cvs2;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent(
    "cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent(
    "cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // we need dg to use correct scaling of basis functions
  Teuchos::ParameterList dglist =
    plist.sublist("PK operator").sublist("flux operator");
  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dglist, mesh0));

  AnalyticDG04 ana(mesh0, order, true);
  ana.InitialGuess(*dg, p1c, 1.0);

  // create remap object
  MyRemapDG remap(mesh0, mesh1, plist);
  remap.DeformMesh(deform, 1.0);
  remap.InitializeOperators(dg);
  remap.InitializeFaceVelocity();
  remap.InitializeJacobianMatrix();
  remap.InitializeConsistentJacobianDeterminant();

  // initial mass
  double mass0(0.0);
  WhetStone::DenseVector data(nk);
  WhetStone::NumericalIntegration numi(mesh0);

  for (int c = 0; c < ncells_owned; c++) {
    for (int i = 0; i < nk; ++i) data(i) = p1c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
    mass0 += numi.IntegratePolynomialCell(c, poly);
  }
  ana.GlobalOp("sum", &mass0, 1);

  // explicit time integration
  CompositeVector p1aux(*p1);
  Explicit_TI::RK<CompositeVector> rk(remap, rk_method, p1aux);

  remap.NonConservativeToConservative(0.0, *p1, p1aux);

  double t(0.0), tend(1.0);
  while (t < tend - dt / 2) {
    // remap.ApplyLimiter(t, p1aux);
    rk.TimeStep(t, dt, p1aux, *p1);
    *p1aux.ViewComponent("cell") = *p1->ViewComponent("cell");

    t += dt;
    remap.CollectStatistics(t, *p1);
  }

  remap.ConservativeToNonConservative(1.0, *p1, p2);

  // calculate error in the new basis
  std::vector<int> dirs;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  CompositeVectorSpace cvs3;
  cvs3.SetMesh(mesh1)->SetGhosted(true)->AddComponent(
    "cell", AmanziMesh::CELL, 1);

  CompositeVector q2(p2);
  Epetra_MultiVector& q2c = *q2.ViewComponent("cell");
  q2c = p2c;

  double pnorm, l2_err, inf_err, l20_err, inf0_err;
  ana.ComputeCellErrorRemap(
    *dg, p2c, tend, 0, mesh1, pnorm, l2_err, inf_err, l20_err, inf0_err);

  CHECK(l2_err < 0.12 / (order + 1));

  if (getRank == 0) {
    printf("nx=%3d (orig) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n",
           nx,
           l20_err,
           l2_err,
           inf0_err,
           inf_err);
  }

  // optional projection on the space of polynomials
  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc0 = mesh0->cell_centroid(c);
    const AmanziGeometry::Point& xc1 = mesh1->cell_centroid(c);

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

    if (order > 0 && order < 3 && dim == 2) {
      poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
      remap.maps()->ProjectPolynomial(c, poly);
      poly.ChangeOrigin(mesh1->cell_centroid(c));
      for (int i = 0; i < nk; ++i) q2c[i][c] = poly(i);
    }
  }

  ana.ComputeCellErrorRemap(
    *dg, q2c, tend, 1, mesh1, pnorm, l2_err, inf_err, l20_err, inf0_err);

  if (getRank == 0) {
    printf("nx=%3d (proj) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n",
           nx,
           l20_err,
           l2_err,
           inf0_err,
           inf_err);
  }

  // concervation errors: mass and volume (CGL)
  auto& det = remap.det();
  double area(0.0), area1(0.0), mass1(0.0), gcl_err(0.0), gcl_inf(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, det[c][0]);
    double vol2 = mesh1->cell_volume(c, false);

    area += vol1;
    area1 += vol2;

    double err = std::fabs(vol1 - vol2);
    gcl_inf = std::max(gcl_inf, err / vol1);
    gcl_err += err;

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

    int quad_order = det[c][0].order() + poly.order();

    if (map_name == "PEM") {
      AmanziMesh::Entity_ID_List faces, nodes;
      mesh0->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      std::vector<AmanziGeometry::Point> xy(3);
      xy[0] = mesh0->cell_centroid(c);

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        mesh0->face_get_nodes(f, &nodes);
        mesh0->node_get_coordinates(nodes[0], &(xy[1]));
        mesh0->node_get_coordinates(nodes[1], &(xy[2]));

        std::vector<const WhetStone::WhetStoneFunction*> polys(2);
        polys[0] = &det[c][n];
        polys[1] = &poly;
        mass1 += numi.IntegrateFunctionsSimplex(xy, polys, quad_order);
      }
    } else {
      WhetStone::Polynomial tmp(det[c][0]);
      tmp.ChangeOrigin(mesh0->cell_centroid(c));
      poly *= tmp;
      mass1 += numi.IntegratePolynomialCell(c, poly);
    }
  }

  // parallel collective operations
  ana.GlobalOp("sum", &area, 1);
  ana.GlobalOp("sum", &area1, 1);
  ana.GlobalOp("sum", &mass1, 1);
  ana.GlobalOp("sum", &gcl_err, 1);
  ana.GlobalOp("max", &gcl_inf, 1);

  if (getRank == 0) {
    printf("Conservation: dMass=%10.4g  dVol=%10.6g  dVolLinear=%10.6g\n",
           mass1 - mass0,
           1.0 - area,
           1.0 - area1);
    printf("GCL: L1=%12.8g  Inf=%12.8g\n", gcl_err, gcl_inf);
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh1, true, false);

  io.InitializeCycle(t, 1);
  io.WriteVector(*p2c(0), "solution");
  io.WriteVector(*q2c(0), "solution-prj");
  io.FinalizeCycle();
}

TEST(REMAP_DUAL_2D)
{
  double dT(0.1);
  auto rk_method = Amanzi::Explicit_TI::heun_euler;
  RemapTestsDualRK(rk_method, "FEM", "", 10, 10, 0, dT);

  RemapTestsDualRK(rk_method, "VEM", "test/median15x16.exo", 0, 0, 0, dT / 2);
  // RemapTestsDualRK(rk_method, "VEM", "", 5,5,5, dT);

  /*
  double dT(0.01);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 1;
  RemapTestsDualRK(rk_method, maps, "",  16, 16,0, dT,    deform);
  RemapTestsDualRK(rk_method, maps, "",  32, 32,0, dT/2,  deform);
  RemapTestsDualRK(rk_method, maps, "",  64, 64,0, dT/4,  deform);
  RemapTestsDualRK(rk_method, maps, "", 128,128,0, dT/8,  deform);
  RemapTestsDualRK(rk_method, maps, "", 256,256,0, dT/16, deform);
  */

  /*
  double dT(0.01);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 4;
  RemapTestsDualRK(rk_method, maps, "test/median15x16.exo",    16,0,0, dT,
  deform); RemapTestsDualRK(rk_method, maps, "test/median32x33.exo",    32,0,0,
  dT/2, deform); RemapTestsDualRK(rk_method, maps, "test/median63x64.exo",
  64,0,0, dT/4, deform); RemapTestsDualRK(rk_method, maps,
  "test/median127x128.exo", 128,0,0, dT/8, deform); RemapTestsDualRK(rk_method,
  maps, "test/median255x256.exo", 256,0,0, dT/16,deform);
  */

  /*
  double dT(0.05);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 4;
  RemapTestsDualRK(rk_method, maps, "test/mesh_poly20x20.exo",    20,0,0, dT,
  deform); RemapTestsDualRK(rk_method, maps, "test/mesh_poly40x40.exo", 40,0,0,
  dT/2, deform); RemapTestsDualRK(rk_method, maps, "test/mesh_poly80x80.exo",
  80,0,0, dT/4, deform); RemapTestsDualRK(rk_method, maps,
  "test/mesh_poly160x160.exo", 160,0,0, dT/8, deform);
  */

  /*
  double dT(0.05);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 4;
  RemapTestsDualRK(rk_method, maps, "test/random10.exo", 10,0,0, dT,   deform);
  RemapTestsDualRK(rk_method, maps, "test/random20.exo", 20,0,0, dT/2, deform);
  RemapTestsDualRK(rk_method, maps, "test/random40.exo", 40,0,0, dT/4, deform);
  */

  /*
  double dT(0.025);
  auto rk_method = rk_method;
  std::string maps = "PEM";
  RemapTestsDualRK(rk_method, maps, "test/triangular8.exo",  0,0,0, dT);
  RemapTestsDualRK(rk_method, maps, "test/triangular16.exo", 0,0,0, dT/2);
  RemapTestsDualRK(rk_method, maps, "test/triangular32.exo", 0,0,0, dT/4);
  RemapTestsDualRK(rk_method, maps, "test/triangular64.exo", 0,0,0, dT/8);
  RemapTestsDualRK(rk_method, maps, "test/triangular128.exo",0,0,0, dT/16);
  */
}
