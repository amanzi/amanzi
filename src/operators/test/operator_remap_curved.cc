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
#include "MeshCurved.hh"
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
            Teuchos::ParameterList& plist, double T1)
    : RemapDG_Tests<AnalyticDG04>(mesh0, mesh1, plist), T1_(T1), tini_(0.0){};
  ~MyRemapDG(){};

  // create basic structures at time zero
  void Init(const Teuchos::RCP<WhetStone::DG_Modal> dg);
  void ReInit(double tini);

  // geometric tools
  virtual void
  DynamicJacobianMatrix(int c, double t, const WhetStone::MatrixPolynomial& J,
                        WhetStone::MatrixPolynomial& Jt) override;
  virtual void DynamicFaceVelocity(double t) override;
  virtual void DynamicCellVelocity(double t) override;

  // mesh deformation from time 0 to t
  virtual void DeformMesh(int deform, double t) override;

  // miscalleneous other tools
  virtual double global_time(double t) override { return tini_ + t * T1_; }

  // access
  const std::vector<WhetStone::VectorPolynomial> det() const { return *det_; }
  const std::shared_ptr<WhetStone::MeshMaps> maps() const { return maps_; }

 private:
  double T1_, tini_;
  std::vector<WhetStone::MatrixPolynomial> J0_;
  std::vector<WhetStone::VectorPolynomial> velf_vec0_;
};


/* *****************************************************************
 * Initialization of remap.
 ***************************************************************** */
void
MyRemapDG::Init(const Teuchos::RCP<WhetStone::DG_Modal> dg)
{
  InitializeOperators(dg);
  InitializeFaceVelocity();
  InitializeJacobianMatrix();

  velf_vec0_.resize(nfaces_wghost_);
  for (int f = 0; f < nfaces_wghost_; ++f) {
    velf_vec0_[f].Reshape(dim_, dim_, 0, true);
  }

  J0_.resize(ncells_owned_);
  for (int c = 0; c < ncells_owned_; ++c) {
    J0_[c].Reshape(dim_, dim_, dim_, 0, true);
    J0_[c].set_origin(mesh0_->cell_centroid(c));
  }
}


/* *****************************************************************
 * TBW
 ***************************************************************** */
void
MyRemapDG::ReInit(double tini)
{
  for (int f = 0; f < nfaces_wghost_; ++f) velf_vec0_[f] += velf_vec_[f];

  for (int c = 0; c < ncells_owned_; ++c) J0_[c] += J_[c];

  InitializeOperators(dg_);
  InitializeFaceVelocity();

  // adjust new velocities for interval [tini, tend]
  for (int f = 0; f < nfaces_wghost_; ++f) velf_vec_[f] -= velf_vec0_[f];

  InitializeJacobianMatrix();

  tini_ = tini;
}


/* *****************************************************************
 * Calculates various geometric quantaties on intermediate meshes.
 ***************************************************************** */
void
MyRemapDG::DynamicJacobianMatrix(int c, double t,
                                 const WhetStone::MatrixPolynomial& J,
                                 WhetStone::MatrixPolynomial& Jt)
{
  Jt = J0_[c] + t * J;

  for (int i = 0; i < dim_; ++i) Jt(i, i)(0) += 1.0;
}


/* *****************************************************************
 * Calculate face co-velocity in reference coordinates
 ***************************************************************** */
void
MyRemapDG::DynamicFaceVelocity(double t)
{
  WhetStone::VectorPolynomial cn; // cn = j J^{-t} N dA

  for (int f = 0; f < nfaces_wghost_; ++f) {
    WhetStone::VectorPolynomial tmp = velf_vec0_[f] + t * velf_vec_[f];
    maps_->NansonFormula(f, tmp, cn);
    (*velf_)[f] = velf_vec_[f] * cn;
  }
}


/* *****************************************************************
 * Cell co-velocity in reference coordinates and Jacobian determinant
 ***************************************************************** */
void
MyRemapDG::DynamicCellVelocity(double t)
{
  WhetStone::MatrixPolynomial Jt, C;
  for (int c = 0; c < ncells_owned_; ++c) {
    DynamicJacobianMatrix(c, t, J_[c], Jt);
    maps_->Determinant(Jt, (*det_)[c]);
    maps_->Cofactors(Jt, C);

    // cell-based pseudo velocity -C^t u
    C.elementWiseMultiply(uc_[c], (*velc_)[c], true);
    (*velc_)[c] *= -1.0;
  }
}


/* *****************************************************************
 * Deform mesh1
 ***************************************************************** */
void
MyRemapDG::DeformMesh(int deform, double t)
{
  RemapDG_Tests<AnalyticDG04>::DeformMesh(deform, t);

  if (order_ > 1) {
    int nfaces =
      mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    auto ho_nodes0 =
      std::make_shared<std::vector<AmanziGeometry::Point_List>>(nfaces);
    auto ho_nodes1 =
      std::make_shared<std::vector<AmanziGeometry::Point_List>>(nfaces);

    for (int f = 0; f < nfaces; ++f) {
      const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
      (*ho_nodes0)[f].push_back(xf);
      (*ho_nodes1)[f].push_back(DeformNode(deform, t, xf));
    }
    auto tmp = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh0_);
    Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(tmp)->set_face_ho_nodes(
      ho_nodes0);
    Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(mesh1_)->set_face_ho_nodes(
      ho_nodes1);
  }
}

} // namespace Amanzi


/* *****************************************************************
 * Remap of polynomilas in two dimensions. Explicit time scheme.
 * Dual formulation places gradient and jumps on a test function.
 ***************************************************************** */
void
RemapTestsCurved(const Amanzi::Explicit_TI::method_t& rk_method,
                 std::string map_name, std::string file_name, int nx, int ny,
                 int nz, double dt0, int deform = 1, int nloop = 1,
                 double T1 = 1.0)
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

  const auto& flux_list = plist.sublist("PK operator").sublist("flux operator");
  int order = flux_list.get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", map_name);

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  const auto& limiter_list = plist.sublist("limiter");
  int vel_order = map_list.get<int>("method order");

  if (getRank == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string limiter = limiter_list.get<std::string>("limiter");
    std::string stencil = limiter_list.get<std::string>("limiter stencil");

    std::cout << "\nTest: " << dim << "D remap:"
              << " mesh=" << ((ny == 0) ? file_name : "square")
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order
              << ", map=" << map_list.get<std::string>("map name")
              << ", flux=" << flux_list.get<std::string>("flux formula")
              << std::endl;

    std::cout << "      map: order=" << vel_order
              << ", projector=" << vel_projector << ", method=\"" << vel_method
              << "\"" << std::endl;

    std::cout << "      limiter: " << limiter << ", stencil=\"" << stencil
              << "\"" << std::endl;
  }

  // create initial mesh
  auto mlist = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("mesh")));
  Teuchos::RCP<MeshCurved> mesh0, mesh1;

  if (dim == 2 && ny != 0) {
    mesh0 =
      Teuchos::rcp(new MeshCurved(0.0, 0.0, 1.0, 1.0, nx, ny, comm, mlist));
    mesh1 =
      Teuchos::rcp(new MeshCurved(0.0, 0.0, 1.0, 1.0, nx, ny, comm, mlist));
  } else if (ny == 0) {
    mesh0 = Teuchos::rcp(new MeshCurved(file_name, comm, mlist));
    mesh1 = Teuchos::rcp(new MeshCurved(file_name, comm, mlist));
  }

  int ncells_owned =
    mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost =
    mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create and initialize cell-based field
  CompositeVectorSpace cvs1, cvs2;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent(
    "cell", AmanziMesh::CELL, nk);
  auto p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent(
    "cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // we need dg to use correct scaling of basis functions
  Teuchos::ParameterList dg_list =
    plist.sublist("PK operator").sublist("flux operator");
  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dg_list, mesh0));

  AnalyticDG04 ana(mesh0, order, true);
  // ana.set_shapes(true, true, false);
  ana.InitialGuess(*dg, p1c, 1.0);

  // visualize initial solution
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh1, true, false));

  p2c = *p1->ViewComponent("cell");

  io->InitializeCycle(0.0, 0);
  io->WriteVector(*p2c(0), "solution");
  io->FinalizeCycle();

  // initial mass
  double mass0(0.0);
  WhetStone::NumericalIntegration numi(mesh0);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p1c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
    mass0 += numi.IntegratePolynomialCell(c, poly);
  }
  ana.GlobalOp("sum", &mass0, 1);

  // create remap object
  MyRemapDG remap(mesh0, mesh1, plist, T1);
  remap.DeformMesh(deform, T1);
  remap.Init(dg);
  remap.set_dt_output(0.1);

  // work in progress on boundary conditions for remap object
  std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0);

  // explicit time integration
  CompositeVector p1aux(*p1);
  Explicit_TI::RK<CompositeVector> rk(remap, rk_method, p1aux);

  remap.NonConservativeToConservative(0.0, *p1, p1aux);

  double dt, t(0.0), tend(0.0);
  for (int iloop = 0; iloop < nloop; ++iloop) {
    dt = dt0;
    tend += T1;
    if (iloop > 0) {
      remap.DeformMesh(deform, tend);
      remap.ReInit(tend - T1);
      t = 0.0;
    }

    // run iterations on pseudo-time interval (0.0, 1.0)
    while (t < 1.0 - dt / 2) {
      rk.TimeStep(t, dt, p1aux, *p1);
      *p1aux.ViewComponent("cell") = *p1->ViewComponent("cell");

      t += dt;
      dt = std::min(dt0, 1.0 - t);
      remap.CollectStatistics(t, *p1);
    }

    remap.ConservativeToNonConservative(t, *p1, p2);

    // visualize solution on mesh1
    // io = Teuchos::rcp(new OutputXDMF(iolist, mesh1, true, false));
    io->InitializeCycle(t, iloop + 1);
    io->WriteVector(*p2c(0), "solution");
    if (order > 0) {
      io->WriteVector(*p2c(1), "gradx");
      io->WriteVector(*p2c(2), "grady");
    }
    if (order > 1) {
      io->WriteVector(*p2c(3), "hesxx");
      io->WriteVector(*p2c(4), "hesxy");
      io->WriteVector(*p2c(5), "hesyy");
    }
    io->FinalizeCycle();
  }

  // calculate error in the new basis
  Entity_ID_List nodes;
  std::vector<int> dirs;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  CompositeVectorSpace cvs3;
  cvs3.SetMesh(mesh1)->SetGhosted(true)->AddComponent(
    "cell", AmanziMesh::CELL, 1);

  double pnorm, l2_err, inf_err, l20_err, inf0_err;
  ana.ComputeCellErrorRemap(
    *dg, p2c, tend, 0, mesh1, pnorm, l2_err, inf_err, l20_err, inf0_err);

  CHECK(l2_err < 0.2 / (order + 1));

  if (getRank == 0) {
    printf("nx=%3d (orig) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n",
           nx,
           l20_err,
           l2_err,
           inf0_err,
           inf_err);
  }

  // optional projection on the space of polynomials
  CompositeVector q2(p2);
  Epetra_MultiVector& q2c = *q2.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc0 = mesh0->cell_centroid(c);
    const AmanziGeometry::Point& xc1 = mesh1->cell_centroid(c);

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];

    if (order > 0 && order < 3 && dim == 2) {
      auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
      remap.maps()->ProjectPolynomial(c, poly);
      poly.ChangeOrigin(mesh1->cell_centroid(c));
      for (int i = 0; i < nk; ++i) q2c[i][c] = poly(i);
    }
  }

  // conservation errors: mass and volume (CGL)
  double area(0.0), area0(0.0), area1(0.0);
  double mass1(0.0), gcl_err(0.0), gcl_inf(0.0);
  auto& det = remap.det();

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, det[c][0]);
    double vol2 = mesh1->cell_volume(c, false);

    area += vol1;
    area0 += mesh0->cell_volume_linear(c);
    area1 += mesh1->cell_volume_linear(c);

    double err = std::fabs(vol1 - vol2);
    gcl_inf = std::max(gcl_inf, err / vol1);
    gcl_err += err;

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

    int quad_order = det[c][0].order() + poly.order();

    WhetStone::Polynomial tmp(det[c][0]);
    tmp.ChangeOrigin(mesh0->cell_centroid(c));
    poly *= tmp;
    mass1 += numi.IntegratePolynomialCell(c, poly);
  }

  // parallel collective operations
  ana.GlobalOp("sum", &area, 1);
  ana.GlobalOp("sum", &area0, 1);
  ana.GlobalOp("sum", &area1, 1);
  ana.GlobalOp("sum", &mass1, 1);
  ana.GlobalOp("sum", &gcl_err, 1);
  ana.GlobalOp("max", &gcl_inf, 1);

  if (getRank == 0) {
    printf("Conservation: dMass=%10.4g  dVolume=%10.6g  dVolLinear=%10.6g\n",
           mass1 - mass0,
           area1 - area,
           area0 - area1);
    printf("GCL: L1=%12.8g  Inf=%12.8g\n", gcl_err, gcl_inf);
  }
}

TEST(REMAP_CURVED_2D)
{
  int nloop = 2;
  double dT(0.1), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::heun_euler;
  std::string maps = "VEM";
  int deform = 1;
  RemapTestsCurved(rk_method, maps, "", 8, 8, 0, dT, deform, nloop, T1);
  // RemapTestsCurved(rk_method, maps, "test/circle_quad10.exo", 10,0,0, 0.1, 6,
  // 40, 0.025);

  /*
  int nloop = 40;
  double dT(0.0025 * nloop), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 6;
  RemapTestsCurved(rk_method, maps, "test/circle_quad10.exo", 10,0,0, dT,
  deform, nloop, T1); RemapTestsCurved(rk_method, maps,
  "test/circle_quad20.exo", 20,0,0, dT/2, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "test/circle_poly40.exo", 40,0,0, dT/4,
  deform, nloop, T1); RemapTestsCurved(rk_method, maps,
  "test/circle_poly80.exo", 80,0,0, dT/8, deform, nloop, T1);
  */

  /*
  int nloop = 5;
  double dT(0.02 * nloop), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 1;
  RemapTestsCurved(rk_method, maps, "",  16, 16,0, dT,   deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "",  32, 32,0, dT/2, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "",  64, 64,0, dT/4, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "", 128,128,0, dT/8, deform, nloop, T1);
  */

  /*
  int nloop = 1;
  double dT(0.01 * nloop), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 4;
  RemapTestsCurved(rk_method, maps, "test/median15x16.exo",    16,0,0, dT,
  deform, nloop, T1); RemapTestsCurved(rk_method, maps, "test/median32x33.exo",
  32,0,0, dT/2, deform, nloop, T1); RemapTestsCurved(rk_method, maps,
  "test/median63x64.exo",    64,0,0, dT/4, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "test/median127x128.exo", 128,0,0, dT/8,
  deform, nloop, T1);
  */
}
