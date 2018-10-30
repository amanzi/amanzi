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
#include "RemapDG.hh"

#include "AnalyticDG04.hh"

namespace Amanzi {

class MyRemapDG : public RemapDG<AnalyticDG04> {
 public:
  MyRemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
            const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
            Teuchos::ParameterList& plist, double T1)
    : RemapDG<AnalyticDG04>(mesh0, mesh1, plist),
      T1_(T1),
      tini_(0.0),
      tl2_(0.0) {};
  ~MyRemapDG() {};

  // create basic structures at time zero
  virtual void Init() override;
  void ReInit(double tini);

  // geometric tools
  virtual void Jacobian(int c, double t, const WhetStone::MatrixPolynomial& J,
                        WhetStone::MatrixPolynomial& Jt) override;
  virtual void UpdateGeometricQuantities(double t) override;

  // mesh deformation from time 0 to t
  virtual void DeformMesh(int deform, double t) override;

  // miscalleneous other tools
  void ChangeVariables(double t, const CompositeVector& p1, CompositeVector& p2, bool flag);
  double L2Norm(double t, double dt, const CompositeVector& p1);
  virtual double global_time(double t) override { return tini_ + t * T1_; }

  // access 
  const std::vector<WhetStone::VectorPolynomial> jac() const { return *jac_; }
  const std::shared_ptr<WhetStone::MeshMaps> maps() const { return maps_; }

 private:
  double T1_, tini_, tl2_;
  std::vector<WhetStone::MatrixPolynomial> J0_;
  std::vector<WhetStone::VectorPolynomial> velf_vec0_;
};


/* *****************************************************************
* Initialization of remap.
***************************************************************** */
void MyRemapDG::Init()
{
  InitPrimary();
  InitSecondary();

  velf_vec0_.resize(nfaces_wghost_);
  for (int f = 0; f < nfaces_wghost_; ++f) {
    velf_vec0_[f].Reshape(dim_, dim_, 0, true);
  }

  J0_.resize(ncells_wghost_);
  for (int c = 0; c < ncells_wghost_; ++c) {
    J0_[c].Reshape(dim_, dim_, dim_, 0, true);
    J0_[c].set_origin(mesh0_->cell_centroid(c));
  }
}


/* *****************************************************************
* TBW
***************************************************************** */
void MyRemapDG::ReInit(double tini)
{
  for (int f = 0; f < nfaces_wghost_; ++f)
    velf_vec0_[f] += velf_vec_[f];

  for (int c = 0; c < ncells_wghost_; ++c)
    J0_[c] += J_[c];

  InitPrimary();

  // adjust new velocities for interval [tini, tend]
  for (int f = 0; f < nfaces_wghost_; ++f)
    velf_vec_[f] -= velf_vec0_[f];

  InitSecondary();

  tini_ = tini;
  tl2_ = 0.0;
}


/* *****************************************************************
* Calculates various geometric quantaties on intermediate meshes.
***************************************************************** */
void MyRemapDG::Jacobian(
    int c, double t, const WhetStone::MatrixPolynomial& J, WhetStone::MatrixPolynomial& Jt)
{
  Jt = J0_[c] + t * J;

  for (int i = 0; i < dim_; ++i)
    Jt(i, i)(0) += 1.0;
}


/* *****************************************************************
* Calculates various geometric quantaties on intermediate meshes.
***************************************************************** */
void MyRemapDG::UpdateGeometricQuantities(double t)
{
  // cn = j J^{-t} N dA
  WhetStone::VectorPolynomial cn;

  for (int f = 0; f < nfaces_wghost_; ++f) {
    WhetStone::VectorPolynomial tmp = velf_vec0_[f] + t * velf_vec_[f]; 
    maps_->NansonFormula(f, tmp, cn);
    (*velf_)[f] = velf_vec_[f] * cn;
  }

  WhetStone::MatrixPolynomial Jt, C;
  for (int c = 0; c < ncells_wghost_; ++c) {
    Jacobian(c, t, J_[c], Jt);
    maps_->Determinant(Jt, (*jac_)[c]);
    maps_->Cofactors(Jt, C);
    
    // cell-based pseudo velocity -C^t u 
    C.Multiply(uc_[c], (*velc_)[c], true);
    (*velc_)[c] *= -1.0;
  }
}


/* *****************************************************************
* Deform mesh1
***************************************************************** */
void MyRemapDG::DeformMesh(int deform, double t)
{
  RemapDG<AnalyticDG04>::DeformMesh(deform, t);

  if (order_ == 2) {
    auto ho_nodes0 = std::make_shared<std::vector<AmanziGeometry::Point_List> >(nfaces_wghost_);
    auto ho_nodes1 = std::make_shared<std::vector<AmanziGeometry::Point_List> >(nfaces_wghost_);

    int nfaces = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
    for (int f = 0; f < nfaces; ++f) {
      const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
      (*ho_nodes0)[f].push_back(xf);
      (*ho_nodes1)[f].push_back(DeformNode(deform, t, xf));
    }
    auto tmp = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh0_);
    Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(tmp)->set_face_ho_nodes(ho_nodes0);
    Teuchos::rcp_static_cast<AmanziMesh::MeshCurved>(mesh1_)->set_face_ho_nodes(ho_nodes1);
  }
}


/* *****************************************************************
* Time integration is done in new variables
***************************************************************** */
void MyRemapDG::ChangeVariables(
    double t, const CompositeVector& p1, CompositeVector& p2, bool flag)
{
  UpdateGeometricQuantities(t);
  op_reac_->Setup(jac_);
  op_reac_->UpdateMatrices(Teuchos::null);

  auto global_reac = op_reac_->global_operator();
  if (flag) {
    global_reac->Apply(p1, p2);
  } else {
    auto& matrices = op_reac_->local_matrices()->matrices;
    for (int n = 0; n < matrices.size(); ++n) {
      matrices[n].Inverse();
    }
    global_reac->Apply(p1, p2);
  }
}


/* *****************************************************************
* L2 norm
***************************************************************** */
double MyRemapDG::L2Norm(double t, double dt, const CompositeVector& p1) {
  if (fabs(tl2_ - t) < dt / 2) {
    CompositeVector p2(p1);

    ChangeVariables(t, p1, p2, false);
    p1.Dot(p2, &l2norm_);
    tl2_ += dt_output_;
  }
  return l2norm_;
} 

}  // namespace Amanzi


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapTestsCurved(const Amanzi::Explicit_TI::method_t& rk_method,
                      std::string map_name, std::string file_name,
                      int nx, int ny, int nz, double dt,
                      int deform = 1, int nloop = 1, double T1 = 1.0) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  int dim = (nz == 0) ? 2 : 3;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  // read parameter list
  std::string xmlFileName = "test/operator_remap.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  int order = plist.sublist("PK operator")
                   .sublist("flux operator").get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", map_name);

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  int vel_order = map_list.get<int>("method order");

  if (MyPID == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string map_name = map_list.get<std::string>("map name");
      
    std::cout << "\nTest: " << dim << "D remap, dual formulation:"
              << " mesh=" << ((ny == 0) ? file_name : "square")
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order 
              << ", map=" << map_name << std::endl;

    std::cout << "      map details: order=" << vel_order 
              << ", projector=" << vel_projector 
              << ", method=\"" << vel_method << "\"" << std::endl;
  }

  // create initial mesh
  const auto& partitioner = AmanziMesh::Partitioner_type::ZOLTAN_RCB;
  Teuchos::RCP<MeshCurved> mesh0, mesh1;

  if (dim == 2 && ny != 0) {
    mesh0 = Teuchos::rcp(new MeshCurved(0.0, 0.0, 1.0, 1.0, nx, ny, &comm, partitioner));
    mesh1 = Teuchos::rcp(new MeshCurved(0.0, 0.0, 1.0, 1.0, nx, ny, &comm, partitioner));
  } else if (ny == 0) {
    mesh0 = Teuchos::rcp(new MeshCurved(file_name, &comm, partitioner));
    mesh1 = Teuchos::rcp(new MeshCurved(file_name, &comm, partitioner));
  }

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
  std::string basis = plist.sublist("PK operator")
                           .sublist("flux operator").get<std::string>("dg basis");
  WhetStone::DG_Modal dg(order, mesh0, basis);
  AnalyticDG04 ana(mesh0, order, true);
  ana.InitialGuess(dg, p1c, 1.0);

  // initial mass
  double mass0(0.0);
  WhetStone::NumericalIntegration numi(mesh0);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p1c[i][c];
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order, data);
    mass0 += numi.IntegratePolynomialCell(c, poly);
  }
  ana.GlobalOp("sum", &mass0, 1);

  // create remap object
  MyRemapDG remap(mesh0, mesh1, plist, T1);
  remap.DeformMesh(deform, T1);
  remap.Init();
  remap.set_dt_output((nloop > 3) ? 0.5 : 0.1);

  // explicit time integration
  CompositeVector p1aux(*p1);
  Explicit_TI::RK<CompositeVector> rk(remap, rk_method, p1aux);

  remap.ChangeVariables(0.0, *p1, p1aux, true);

  int nstep(0);
  double t(0.0), tend(0.0);
  for (int iloop = 0; iloop < nloop; ++iloop) {
    tend += T1;
    if (iloop > 0) { 
      remap.DeformMesh(deform, tend);
      remap.ReInit(tend - T1);
      t = 0.0;
    }

    // run iterations on pseudo-time interval (0.0, 1.0)
    while(t < 1.0 - dt/2) {
      remap.L2Norm(t, dt, p1aux);
      rk.TimeStep(t, dt, p1aux, *p1);

      *p1aux.ViewComponent("cell") = *p1->ViewComponent("cell");

      t += dt;
      nstep++;
    }
  }

  remap.ChangeVariables(1.0, *p1, p2, false);

  // calculate error in the new basis
  Entity_ID_List nodes;
  std::vector<int> dirs;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  CompositeVectorSpace cvs3;
  cvs3.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  CompositeVector q2(p2);
  Epetra_MultiVector& q2c = *q2.ViewComponent("cell");
  q2c = p2c;

  double pnorm, l2_err, inf_err, l20_err, inf0_err;
  ana.ComputeCellErrorRemap(dg, p2c, tend, 0, mesh1,
                            pnorm, l2_err, inf_err, l20_err, inf0_err);

  CHECK(l2_err < 0.2 / (order + 1));

  if (MyPID == 0) {
    printf("nx=%3d (orig) L2=%12.8g %12.8g  Inf=%12.8g %12.8g\n", 
        nx, l20_err, l2_err, inf0_err, inf_err);
  }

  // conservation errors: mass and volume (CGL)
  double area(0.0), area1(0.0), mass1(0.0), gcl_err(0.0), gcl_inf(0.0);
  auto& jac = remap.jac();

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, jac[c][0]);
    double vol2 = mesh1->cell_volume(c);

    area += vol1;
    area1 += vol2;

    double err = std::fabs(vol1 - vol2);
    gcl_inf = std::max(gcl_inf, err / vol1);
    gcl_err += err;

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];
    auto poly = dg.cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

    int quad_order = jac[c][0].order() + poly.order();

    WhetStone::Polynomial tmp(jac[c][0]);
    tmp.ChangeOrigin(mesh0->cell_centroid(c));
    poly *= tmp;
    mass1 += numi.IntegratePolynomialCell(c, poly);
  }

  // parallel collective operations
  ana.GlobalOp("sum", &area, 1);
  ana.GlobalOp("sum", &area1, 1);
  ana.GlobalOp("sum", &mass1, 1);
  ana.GlobalOp("sum", &gcl_err, 1);
  ana.GlobalOp("max", &gcl_inf, 1);

  if (MyPID == 0) {
    printf("Conservation: dMass=%10.4g  dVolume=%10.6g  dVolLinear=%10.6g\n",
           mass1 - mass0, 1.0 - area, 1.0 - area1);
    printf("GCL: L1=%12.8g  Inf=%12.8g\n", gcl_err, gcl_inf);
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh1, true, false);

  io.InitializeCycle(t, nstep);
  io.WriteVector(*p2c(0), "remapped");
  io.WriteVector(*q2c(0), "remapped-prj");
  io.FinalizeCycle();
}

TEST(REMAP_CURVED_2D) {
  int nloop = 2;
  double dT(0.1), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::heun_euler;
  std::string maps = "VEM";
  int deform = 1;
  RemapTestsCurved(rk_method, maps, "", 8,8,0, dT, deform, nloop, T1);

  /*
  int nloop = 5;
  double dT(0.01 * nloop), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 4;
  RemapTestsCurved(rk_method, maps, "",  16, 16,0, dT,   deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "",  32, 32,0, dT/2, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "",  64, 64,0, dT/4, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "", 128,128,0, dT/8, deform, nloop, T1);
  */

  /*
  int nloop = 5;
  double dT(0.01 * nloop), T1(1.0 / nloop);
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;
  std::string maps = "VEM";
  int deform = 4;
  RemapTestsCurved(rk_method, maps, "test/median15x16.exo",    16,0,0, dT,   deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "test/median32x33.exo",    32,0,0, dT/2, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "test/median63x64.exo",    64,0,0, dT/4, deform, nloop, T1);
  RemapTestsCurved(rk_method, maps, "test/median127x128.exo", 128,0,0, dT/8, deform, nloop, T1);
  */
}

