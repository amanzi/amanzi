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
#include "TreeVector.hh"

// Amanzi::Operators
#include "MeshDeformation.hh"
#include "RemapDG.hh"

#include "AnalyticDG04.hh"

namespace Amanzi {

class MyRemapDG : public Operators::RemapDG<TreeVector> {
 public:
  MyRemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
            const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
            Teuchos::ParameterList& plist)
    : RemapDG<TreeVector>(mesh0, mesh1, plist),
      tprint_(0.0),
      dt_output_(0.1),
      l2norm_(-1.0),
      T1_(1.0),
      tini_(0.0) {};
  ~MyRemapDG() {};

  // time control
  double global_time(double t) { return tini_ + t * T1_; }

  // tools
  // -- mass on mesh0
  double InitialMass(const TreeVector& p1, int order);
  // -- statistics
  void CollectStatistics(double t, const TreeVector& u);

  // access 
  const std::vector<WhetStone::Polynomial> jac() const { return *jac_; }

 public:
  double tprint_, dt_output_, l2norm_;

 private:
  double T1_, tini_;
};


/* *****************************************************************
* Compute initial mass: partial specialization
***************************************************************** */
double MyRemapDG::InitialMass(const TreeVector& p1, int order)
{
  const Epetra_MultiVector& p1c = *p1.SubVector(0)->Data()->ViewComponent("cell", false);
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
void MyRemapDG::CollectStatistics(double t, const TreeVector& u)
{
  double tglob = global_time(t);
  if (tglob >= tprint_) {
    op_reac_->Setup(det_, false);
    op_reac_->UpdateMatrices(t);
    auto& matrices = op_reac_->local_op()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

    auto& rhs = *op_reac_->global_operator()->rhs();
    op_reac_->global_operator()->Apply(*u.SubVector(0)->Data(), rhs);
    rhs.Dot(*u.SubVector(0)->Data(), &l2norm_);

    Epetra_MultiVector& xc = *rhs.ViewComponent("cell");
    int nk = xc.NumVectors();
    double xmax[nk], xmin[nk];
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (mesh0_->get_comm()->MyPID() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  umax/umin: %9.5g %9.5g\n",
             tglob, l2norm_, nfun_, sharp_, xmax[0], xmin[0]);
    }

    tprint_ += dt_output_;
    sharp_ = 0.0;
  } 
}

}  // namespace Amanzi


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapGCL(const Amanzi::Explicit_TI::method_t& rk_method,
              std::string file_name,
              int nx, int ny, int nz, double dt, int deform = 1) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  int dim = (nz == 0) ? 2 : 3;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // read parameter list
  std::string xmlFileName = "test/operator_remap.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  int order = plist.sublist("PK operator")
                   .sublist("flux operator")
                   .sublist("schema").get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", "VEM");

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  int vel_order = map_list.get<int>("method order");

  if (MyPID == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
      
    std::cout << "\nTest: " << dim << "D remap:"
              << " mesh=" << ((file_name == "") ? "structured" : file_name)
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order << std::endl;

    std::cout << "      map details: order=" << vel_order 
              << ", projector=" << vel_projector 
              << ", method=\"" << vel_method << "\"" << std::endl;

    std::cout << "      RK method: " << (int)rk_method << std::endl;
  }

  // create two meshes
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

  // create and initialize cell-based fields
  auto cvs1 = Teuchos::rcp(new CompositeVectorSpace());
  cvs1->SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);

  Teuchos::RCP<TreeVectorSpace> tvs0 = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> tvs1 = Teuchos::rcp(new TreeVectorSpace());
  tvs0->SetData(cvs1);
  tvs1->PushBack(tvs0);
  tvs1->PushBack(tvs0);

  Teuchos::RCP<TreeVector> p1 = Teuchos::rcp(new TreeVector(*tvs1));
  Epetra_MultiVector& p1c = *p1->SubVector(0)->Data()->ViewComponent("cell", true);
  Epetra_MultiVector& j1c = *p1->SubVector(1)->Data()->ViewComponent("cell", true);

  auto cvs2 = Teuchos::rcp(new CompositeVectorSpace());
  cvs2->SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);

  Teuchos::RCP<TreeVectorSpace> tvs2 = Teuchos::rcp(new TreeVectorSpace());
  tvs0->SetData(cvs2);
  tvs2->PushBack(tvs0);
  tvs2->PushBack(tvs0);

  Teuchos::RCP<TreeVector> p2 = Teuchos::rcp(new TreeVector(*tvs2));
  Epetra_MultiVector& p2c = *p2->SubVector(0)->Data()->ViewComponent("cell");

  // we need dg to use correct scaling of basis functions
  Teuchos::ParameterList dglist = plist.sublist("PK operator")
                                       .sublist("flux operator").sublist("schema");
  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dglist, mesh0));

  AnalyticDG04 ana(mesh0, order, true);
  ana.InitialGuess(*dg, p1c, 1.0);
  j1c.PutScalar(0.0);
  j1c(0)->PutScalar(1.0);

  // create remap object
  MyRemapDG remap(mesh0, mesh1, plist);
  DeformMeshCurved(mesh1, deform, 1.0, mesh0, order);
  remap.InitializeOperators(dg);
  if (MyPID == 0) std::cout << "Computing static data on mesh scheleton...\n";
  remap.StaticEdgeFaceVelocities();
  remap.StaticFaceCoVelocity();
  if (MyPID == 0) std::cout << "Computing static data in mesh cells...\n";
  remap.StaticCellVelocity();
  remap.StaticCellCoVelocity();

  // initial mass
  double mass0 = remap.InitialMass(*p1, order);

  // explicit time integration
  Explicit_TI::RK<TreeVector> rk(remap, rk_method, *p1);

  TreeVector p3(*p1);
  auto& sv3 = *p3.SubVector(0)->Data();
  Epetra_MultiVector& p3c = *sv3.ViewComponent("cell", true);
  remap.NonConservativeToConservative(0.0, *p1, p3);

  double t(0.0), tend(1.0);
  while(t < tend - dt/2) {
    // remap.ApplyLimiter(t, p3);
    rk.TimeStep(t, dt, p3, *p1);
    p3 = *p1;

    t += dt;
    remap.CollectStatistics(t, *p1);
  }

  remap.ConservativeToNonConservative(1.0, *p1, *p2);

  // calculate error in the new basis
  std::vector<int> dirs;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  double pnorm, l2_err, inf_err, l20_err, l10_err, inf0_err;
  ana.ComputeCellErrorRemap(*dg, p2c, tend, 0, mesh1,
                            pnorm, l2_err, inf_err, l20_err, l10_err, inf0_err, &p3c);

  CHECK(((dim == 2) ? l2_err : l20_err) < 0.12 / (order + 1));

  if (MyPID == 0) {
    printf("nx=%3d (orig) L1=%12.8g(mean) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n", 
        nx, l10_err, l20_err, l2_err, inf0_err, inf_err);
  }

  // concervation errors: mass and volume (CGL)
  auto& jac = remap.jac();
  double area(0.0), area1(0.0), mass1(0.0), gcl_err(0.0), gcl_inf(0.0);
  WhetStone::NumericalIntegration numi(mesh0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, jac[c]);
    double vol2 = mesh1->cell_volume(c);

    area += vol1;
    area1 += vol2;

    double err = std::fabs(vol1 - vol2);
    gcl_inf = std::max(gcl_inf, err / vol1);
    gcl_err += err;

    WhetStone::DenseVector data(nk);
    for (int i = 0; i < nk; ++i) data(i) = p2c[i][c];
    auto poly = dg->cell_basis(c).CalculatePolynomial(mesh0, c, order, data);

    WhetStone::Polynomial tmp(jac[c]);
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
    printf("Conservation: dMass=%10.4g  dVol=%10.6g  dVolLinear=%10.6g\n",
           mass1 - mass0, 1.0 - area, 1.0 - area1);
    printf("GCL: L1=%12.8g  Inf=%12.8g\n", gcl_err, gcl_inf);
  }
  CHECK_CLOSE(mass0, mass1, 1e-12);

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh1, true, false);

  io.InitializeCycle(t, 1, "");
  io.WriteVector(*p2c(0), "solution", AmanziMesh::CELL);
  io.FinalizeCycle();
}

TEST(REMAP_GEOMETRIC_CONSERVATION_LAW) {
  int deform = 5;
  auto rk_method = Amanzi::Explicit_TI::tvd_3rd_order;

  double dT(0.1);
  RemapGCL(rk_method, "test/median15x16.exo",   16,1,0, dT/2, deform);
  /*
  RemapGCL(rk_method, "test/median32x33.exo",   32,1,0, dT/4, deform);
  RemapGCL(rk_method, "test/median63x64.exo",   64,1,0, dT/8, deform);
  RemapGCL(rk_method, "test/median127x128.exo",128,1,0, dT/16,deform);

  double dT(0.025);
  RemapGCL(rk_method, "test/prism10.exo", 10,1,1, dT,   deform);
  RemapGCL(rk_method, "test/prism20.exo", 20,1,1, dT/2, deform);
  RemapGCL(rk_method, "test/prism40.exo", 40,1,1, dT/4, deform);
  */
}

