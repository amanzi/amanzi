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
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CompositeVector.hh"
#include "DG_Modal.hh"
#include "Explicit_TI_RK.hh"
#include "Mesh.hh"
#include "MeshCurved.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "OutputXDMF.hh"

// Amanzi::Operators
#include "MeshDeformation.hh"
#include "RemapDG.hh"

#include "AnalyticDG04.hh"
#include "AnalyticDG04b.hh"

namespace Amanzi {

class MyRemapDGc : public Operators::RemapDG<CompositeVector> {
 public:
  MyRemapDGc(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
             const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
             Teuchos::ParameterList& plist,
             double T1)
    : RemapDG<CompositeVector>(mesh0, mesh1, plist),
      tprint_(0.0),
      dt_output_(0.1),
      l2norm_(-1.0),
      T1_(T1),
      tini_(0.0){};
  ~MyRemapDGc(){};

  // create basic structures at time zero
  void Init(const Teuchos::RCP<WhetStone::DG_Modal> dg);
  void ReInit(double tini);

  // co-velocities
  virtual void StaticFaceCoVelocity() override;
  virtual void StaticCellCoVelocity() override;

  // time control
  double global_time(double t) { return tini_ + t * T1_; }
  void set_dt_output(double dt) { dt_output_ = dt; }

  // tools
  // -- mass on mesh0
  double InitialMass(const CompositeVector& p1, int order);
  // -- statictics
  void CollectStatistics(double t, const CompositeVector& u);

  // access
  const std::vector<WhetStone::SpaceTimePolynomial> det() const { return *det_; }
  const std::shared_ptr<WhetStone::MeshMaps> maps() const { return maps_; }

 public:
  double tprint_, dt_output_, l2norm_;

 private:
  double T1_, tini_;
  std::vector<WhetStone::VectorPolynomial> velf_vec0_, vele_vec0_;
  std::vector<WhetStone::MatrixPolynomial> J_, J0_;
};


/* *****************************************************************
* Initialization of remap.
***************************************************************** */
void
MyRemapDGc::Init(const Teuchos::RCP<WhetStone::DG_Modal> dg)
{
  if (mesh0_->get_comm()->MyPID() == 0) std::cout << "Computing static data on mesh scheleton...\n";
  InitializeOperators(dg);
  StaticEdgeFaceVelocities();
  StaticCellVelocity();

  velf_vec0_.resize(nfaces_wghost_);
  for (int f = 0; f < nfaces_wghost_; ++f) {
    velf_vec0_[f].Reshape(dim_, dim_, 1, true);
    velf_vec0_[f].set_origin(mesh0_->face_centroid(f));
  }

  if (mesh0_->valid_edges()) {
    vele_vec0_.resize(nedges_wghost_);
    for (int e = 0; e < nedges_wghost_; ++e) {
      vele_vec0_[e].Reshape(dim_, dim_, 1, true);
      vele_vec0_[e].set_origin(mesh0_->edge_centroid(e));
    }
  }

  J0_.resize(ncells_owned_);
  for (int c = 0; c < ncells_owned_; ++c) {
    J0_[c].Reshape(dim_, dim_, dim_, 0, true);
    J0_[c].set_origin(mesh0_->cell_centroid(c));
  }

  if (mesh0_->get_comm()->MyPID() == 0) std::cout << "Computing static data in mesh cells...\n";
  StaticFaceCoVelocity();
  StaticCellCoVelocity();
}


/* *****************************************************************
* TBW
***************************************************************** */
void
MyRemapDGc::ReInit(double tini)
{
  if (mesh0_->get_comm()->MyPID() == 0) std::cout << "Computing static data on mesh scheleton...\n";
  for (int f = 0; f < nfaces_wghost_; ++f) velf_vec0_[f] += velf_vec_[f];

  if (mesh0_->valid_edges()) {
    for (int e = 0; e < nedges_wghost_; ++e) vele_vec0_[e] += vele_vec_[e];
  }

  for (int c = 0; c < ncells_owned_; ++c) J0_[c] += J_[c];

  InitializeOperators(dg_);
  StaticEdgeFaceVelocities();

  // adjust new velocities for interval [tini, tend]
  for (int f = 0; f < nfaces_wghost_; ++f) velf_vec_[f] -= velf_vec0_[f];

  if (mesh0_->valid_edges()) {
    for (int e = 0; e < nedges_wghost_; ++e) vele_vec_[e] -= vele_vec0_[e];
  }

  StaticCellVelocity();

  if (mesh0_->get_comm()->MyPID() == 0) std::cout << "Computing static data in mesh cells...\n";
  StaticFaceCoVelocity();
  StaticCellCoVelocity();

  // re-calculate static matrices
  if (mesh0_->get_comm()->MyPID() == 0) std::cout << "Computing static matrices for operators...\n";
  op_adv_->Setup(velc_, true);
  op_reac_->Setup(det_, true);
  op_flux_->Setup(velf_.ptr(), true);

  tini_ = tini;
}


/* *****************************************************************
* Initialization of space-tim co-velocity v = u * (j J^{-t} N)
***************************************************************** */
void
MyRemapDGc::StaticFaceCoVelocity()
{
  WhetStone::VectorSpaceTimePolynomial cn;
  for (int f = 0; f < nfaces_wghost_; ++f) {
    WhetStone::VectorSpaceTimePolynomial map(dim_, dim_, 1), tmp(dim_, dim_, 0);

    for (int i = 0; i < dim_; ++i) {
      map[i][0] = velf_vec0_[f][i]; // map = u0
      map[i][0](1, i) += 1.0;       // map = x + u0
      map[i][1] = velf_vec_[f][i];  // map = x + u0 + t * (u - u0)

      tmp[i][0] = velf_vec_[f][i];
    }

    maps_->NansonFormula(f, map, cn);
    (*velf_)[f] = tmp * cn;
  }
}


/* *****************************************************************
* Initialization of the constant cell velocity
***************************************************************** */
void
MyRemapDGc::StaticCellCoVelocity()
{
  J_.resize(ncells_owned_);
  for (int c = 0; c < ncells_owned_; ++c) {
    maps_->Jacobian(uc_[c], J_[c]);

    // space-time cell velocity: v = -j J^{-1} u = -C^t u
    WhetStone::MatrixSpaceTimePolynomial Jt(dim_, dim_, dim_, 1), Ct;
    WhetStone::VectorSpaceTimePolynomial tmp(dim_, dim_, 0);

    for (int i = 0; i < dim_; ++i) {
      for (int j = 0; j < dim_; ++j) {
        Jt(i, j)[0] = J0_[c](i, j);
        Jt(i, j)[1] = J_[c](i, j); // Jt = J0 + t * J
      }
      Jt(i, i)[0](0) += 1.0;
      tmp[i][0] = uc_[c][i];
    }

    maps_->Cofactors(Jt, Ct);

    tmp *= -1.0;
    Ct.Multiply(tmp, (*velc_)[c], true);

    maps_->Determinant(Jt, (*det_)[c]);
  }
}


/* *****************************************************************
* Compute initial mass: partial specialization
***************************************************************** */
double
MyRemapDGc::InitialMass(const CompositeVector& p1, int order)
{
  const Epetra_MultiVector& p1c = *p1.ViewComponent("cell", false);
  int nk = p1c.NumVectors();
  int ncells = p1c.MyLength();

  double mass(0.0), mass0;
  WhetStone::DenseVector data(nk);
  WhetStone::NumericalIntegration numi(mesh0_);

  for (int c = 0; c < ncells; c++) {
    for (int i = 0; i < nk; ++i) data(i) = p1c[i][c];
    auto poly = dg_->cell_basis(c).CalculatePolynomial(mesh0_, c, order, data);
    mass += numi.IntegratePolynomialCell(c, poly);
  }

  mesh0_->get_comm()->SumAll(&mass, &mass0, 1);
  return mass0;
}


/* *****************************************************************
* Print statistics using conservative field u
***************************************************************** */
void
MyRemapDGc::CollectStatistics(double t, const CompositeVector& u)
{
  double tglob = global_time(t);
  if (tglob >= tprint_) {
    op_reac_->UpdateMatrices(t);
    auto& matrices = op_reac_->local_op()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();

    auto& rhs = *op_reac_->global_operator()->rhs();
    op_reac_->global_operator()->Apply(u, rhs);
    rhs.Dot(u, &l2norm_);

    Epetra_MultiVector& xc = *rhs.ViewComponent("cell");
    int nk = xc.NumVectors();
    double xmax[nk], xmin[nk], lmax(-1.0), lmin(-1.0), lavg(-1.0);
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (limiter() != Teuchos::null) {
      const auto& lim = *limiter()->limiter();
      lim.MaxValue(&lmax);
      lim.MinValue(&lmin);
      lim.MeanValue(&lavg);
    }

    if (mesh0_->get_comm()->MyPID() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  limiter: %6.3f %6.3f %6.3f  umax/umin: "
             "%9.5g %9.5g\n",
             tglob,
             l2norm_,
             nfun_,
             sharp_,
             lmax,
             lmin,
             lavg,
             xmax[0],
             xmin[0]);
    }

    tprint_ += dt_output_;
    sharp_ = 0.0;
  }
}

} // namespace Amanzi


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
template <class AnalyticDG>
void
RemapTestsCurved(std::string file_name,
                 int nx,
                 int ny,
                 int nz,
                 double dt0,
                 int deform = 1,
                 int nloop = 1,
                 double T1 = 1.0)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  int dim = (nz == 0) ? 2 : 3;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // read and set parameters
  std::string xmlFileName = "test/operator_remap.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  const auto& flux_list = plist.sublist("PK operator").sublist("flux operator");
  int order = flux_list.sublist("schema").get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", "VEM");

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  const auto& limiter_list = plist.sublist("limiter");
  int vel_order = map_list.get<int>("method order");

  if (MyPID == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string limiter = limiter_list.get<std::string>("limiter");
    std::string stencil = limiter_list.get<std::string>("limiter stencil");

    std::cout << "\nTest: " << dim << "D remap curved:"
              << " mesh=" << ((ny == 0) ? file_name : "box mesh") << " deform=" << deform
              << std::endl;

    std::cout << "      discretization: order=" << order
              << ", map=" << map_list.get<std::string>("map name")
              << ", flux=" << flux_list.get<std::string>("flux formula") << std::endl;

    std::cout << "      map: order=" << vel_order << ", projector=" << vel_projector
              << ", method=\"" << vel_method << "\"" << std::endl;

    std::cout << "      limiter: " << limiter << ", stencil=\"" << stencil << "\"" << std::endl;

    std::cout << "      RK method: " << (int)rk_method << std::endl;
  }

  // create initial mesh
  Teuchos::ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(dim, region_list, *comm));

  auto mlist = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("mesh")));
  Teuchos::RCP<MeshCurved> mesh0, mesh1;

  if (file_name != "") {
    bool request_edges = (dim == 3);
    mesh0 = Teuchos::rcp(new MeshCurved(file_name, comm, gm, mlist, true, request_edges));
    mesh1 = Teuchos::rcp(new MeshCurved(file_name, comm, gm, mlist, true, request_edges));
  } else if (dim == 2) {
    mesh0 = Teuchos::rcp(new MeshCurved(0.0, 0.0, 1.0, 1.0, nx, ny, comm, gm, mlist));
    mesh1 = Teuchos::rcp(new MeshCurved(0.0, 0.0, 1.0, 1.0, nx, ny, comm, gm, mlist));
  } else {
    mesh0 = Teuchos::rcp(new MeshCurved(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, comm, gm, mlist));
    mesh1 = Teuchos::rcp(new MeshCurved(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, comm, gm, mlist));
  }
  mesh0->BuildCache();
  mesh1->BuildCache();

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create and initialize cell-based field
  CompositeVectorSpace cvs1, cvs2;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  auto p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // we need dg to use correct scaling of basis functions
  Teuchos::ParameterList dg_list =
    plist.sublist("PK operator").sublist("flux operator").sublist("schema");
  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dg_list, mesh0));

  AnalyticDG ana(mesh0, order, true);
  // ana.set_shapes(true, true, true);
  ana.InitialGuess(*dg, p1c, 1.0);

  // visualize initial solution
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh1, true, false));

  p2c = *p1->ViewComponent("cell");

  io->InitializeCycle(0.0, 0, "");
  io->WriteVector(*p2c(0), "solution", AmanziMesh::CELL);
  io->FinalizeCycle();

  // create remap object
  MyRemapDGc remap(mesh0, mesh1, plist, T1);
  DeformMeshCurved(mesh1, deform, T1, mesh0, order);
  remap.Init(dg);
  remap.set_dt_output(0.1);

  // initial mass
  double mass0 = remap.InitialMass(*p1, order);

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
      DeformMeshCurved(mesh1, deform, tend, mesh0, order);
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
    io->InitializeCycle(t, iloop + 1, "");
    io->WriteVector(*p2c(0), "solution", AmanziMesh::CELL);
    if (order > 0) {
      io->WriteVector(*p2c(1), "gradx", AmanziMesh::CELL);
      io->WriteVector(*p2c(2), "grady", AmanziMesh::CELL);
    }
    if (order > 1) {
      io->WriteVector(*p2c(3), "hesxx", AmanziMesh::CELL);
      io->WriteVector(*p2c(4), "hesxy", AmanziMesh::CELL);
      io->WriteVector(*p2c(5), "hesyy", AmanziMesh::CELL);
    }
    io->FinalizeCycle();
  }

  // calculate error in the new basis
  Entity_ID_List nodes;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  double pnorm, l2_err, inf_err, l20_err, l10_err, inf0_err;
  ana.ComputeCellErrorRemap(
    *dg, p2c, tend, 0, mesh1, pnorm, l2_err, inf_err, l20_err, l10_err, inf0_err);

  CHECK(l20_err < 0.2 / (order + 1));

  if (MyPID == 0) {
    printf("nx=%3d (orig) L1=%12.8g(mean) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n",
           nx,
           l10_err,
           l20_err,
           l2_err,
           inf0_err,
           inf_err);
  }

  // optional projection on the space of polynomials
  CompositeVector q2(p2);
  Epetra_MultiVector& q2c = *q2.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
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
  WhetStone::NumericalIntegration numi(mesh0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, det[c].Value(1.0));
    double vol2 = mesh1->cell_volume(c);

    area += vol1;
    area0 += mesh0->cell_volume_linear(c);
    area1 += mesh1->cell_volume_linear(c);

    double err = std::fabs(vol1 - vol2);
    gcl_inf = std::max(gcl_inf, err / vol1);
    gcl_err += err;

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

    WhetStone::Polynomial tmp(det[c].Value(1.0));
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

  if (MyPID == 0) {
    printf("Conservation: dMass=%10.4g  dVolume=%10.6g  dVolLinear=%10.6g\n",
           mass1 - mass0,
           area1 - area,
           area0 - area1);
    printf("GCL: L1=%12.8g  Inf=%12.8g\n", gcl_err, gcl_inf);
  }
}

TEST(REMAP_CURVED_3D)
{
  int nloop = 1;
  double dT(0.025), T1(1.0 / nloop);
  int deform = 5;
  RemapTestsCurved<AnalyticDG04b>("", 4, 4, 4, dT, deform, nloop, T1);
  // RemapTestsCurved<AnalyticDG04b>("test/prism10.exo", 10,1,1, dT,   deform, nloop, T1);
}

TEST(REMAP_CURVED_2D)
{
  int nloop = 2;
  double dT(0.1), T1(1.0 / nloop);
  int deform = 1;
  RemapTestsCurved<AnalyticDG04>("", 8, 8, 0, dT, deform, nloop, T1);
  // RemapTestsCurved<AnalyticDG04>("test/circle_quad10.exo", 10,0,0, 0.1, 6, 40, 0.025);
}

TEST(REMAP_CURVED_DEV)
{
  /*
  int nloop = 2;
  double dT(0.025), T1(1.0 / nloop);
  int deform = 5;
  RemapTestsCurved<AnalyticDG04b>("test/prism10.exo", 10,1,1, dT,   deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04b>("test/prism20.exo", 20,1,1, dT/2, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04b>("test/prism40.exo", 40,1,1, dT/4, deform, nloop, T1);
  */

  /*
  int nloop = 40;
  double dT(0.0025 * nloop), T1(1.0 / nloop);
  int deform = 6;
  RemapTestsCurved<AnalyticDG04>("test/circle_quad10.exo", 10,0,0, dT,   deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("test/circle_quad20.exo", 20,0,0, dT/2, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("test/circle_poly40.exo", 40,0,0, dT/4, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("test/circle_poly80.exo", 80,0,0, dT/8, deform, nloop, T1);
  */

  /*
  int nloop = 2;
  double dT(0.01 * nloop), T1(1.0 / nloop);
  int deform = 1;
  RemapTestsCurved<AnalyticDG04>("",  16, 16,0, dT,   deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("",  32, 32,0, dT/2, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("",  64, 64,0, dT/4, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("", 128,128,0, dT/8, deform, nloop, T1);
  */

  /*
  int nloop = 1;
  double dT(0.01 * nloop), T1(1.0 / nloop);
  int deform = 5;
  RemapTestsCurved<AnalyticDG04>("test/median15x16.exo",    16,0,0, dT,   deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("test/median32x33.exo",    32,0,0, dT/2, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("test/median63x64.exo",    64,0,0, dT/4, deform, nloop, T1);
  RemapTestsCurved<AnalyticDG04>("test/median127x128.exo", 128,0,0, dT/8, deform, nloop, T1);
  */
}
