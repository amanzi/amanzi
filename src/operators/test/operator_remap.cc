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
#include "CompositeVector.hh"
#include "DG_Modal.hh"
#include "Explicit_TI_RK.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "OutputXDMF.hh"

// Amanzi::Operators
#include "MeshDeformation.hh"
#include "RemapDG.hh"

#include "AnalyticDG04.hh"

namespace Amanzi {

class MyRemapDG : public Operators::RemapDG<CompositeVector> {
 public:
  MyRemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
            const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
            Teuchos::ParameterList& plist)
    : Operators::RemapDG<CompositeVector>(mesh0, mesh1, plist),
      tprint_(0.0),
      l2norm_(-1.0),
      dt_output_(0.1) {};
  ~MyRemapDG() {};

  // time control
  // -- stability condition
  double StabilityCondition();
  // -- time
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
};


/* *****************************************************************
* Rough estimate of the CFL condition.
***************************************************************** */
double MyRemapDG::StabilityCondition()
{
  double dt(1e+99), alpha(0.2), tmp;

  for (int f = 0; f < nfaces_wghost_; ++f) {
    double area = mesh0_->face_area(f);
    const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
    velf_vec_[f].Value(xf).Norm2(&tmp);
    dt = std::min(dt, area / tmp);
  }

  return dt * alpha / (2 * order_ + 1);
}


/* *****************************************************************
* Compute initial mass: partial specialization
***************************************************************** */
double MyRemapDG::InitialMass(const CompositeVector& p1, int order)
{
  const Epetra_MultiVector& p1c = *p1.ViewComponent("cell", false);
  int nk = p1c.NumVectors();
  int ncells = p1c.MyLength();

  double mass(0.0), mass0;
  WhetStone::DenseVector data(nk);
  WhetStone::NumericalIntegration<AmanziMesh::MeshLight> numi(mesh0_);

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
void MyRemapDG::CollectStatistics(double t, const CompositeVector& u)
{
  double tglob = t;
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
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  limiter: %6.3f %6.3f %6.3f  umax/umin: %9.5g %9.5g\n",
             tglob, l2norm_, nfun_, sharp_, lmax, lmin, lavg, xmax[0], xmin[0]);
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
void RemapTestsDualRK(std::string map_name, std::string file_name,
                      int nx, int ny, int nz, double dt,
                      int deform = 1) {
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

  int order = plist.sublist("PK operator")
                   .sublist("flux operator")
                   .sublist("schema").get<int>("method order");

  int nk = WhetStone::PolynomialSpaceDimension(dim, order);

  auto rk_method = Amanzi::Explicit_TI::heun_euler;

  // make modifications to the parameter list
  plist.sublist("maps").set<std::string>("map name", map_name);

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  int vel_order = map_list.get<int>("method order");

  if (MyPID == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string name = map_list.get<std::string>("map name");
      
    std::cout << "\nTest: " << dim << "D remap:"
              << " mesh=" << ((file_name == "") ? "structured" : file_name)
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order 
              << ", map=" << name << std::endl;

    std::cout << "      map details: order=" << vel_order 
              << ", projector=" << vel_projector 
              << ", method=\"" << vel_method << "\"" << std::endl;

    std::cout << "      RK method: " << (int)rk_method << std::endl;
  }

  // create two meshes
  Teuchos::ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(dim, region_list, *comm));

  auto mlist = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("mesh")));
  MeshFactory meshfactory(comm, gm, mlist);
  meshfactory.set_preference(Preference({AmanziMesh::Framework::MSTK}));

  Teuchos::RCP<const Mesh> mesh0;
  Teuchos::RCP<Mesh> mesh1;
  if (file_name != "") {
    bool request_edges = (dim == 3);
    mesh0 = meshfactory.create(file_name, true, request_edges);
    mesh1 = meshfactory.create(file_name, true, request_edges);
  } else if (dim == 2) {
    mesh0 = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, ny);
    mesh1 = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, ny);
  } else { 
    mesh0 = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, true, true);
    mesh1 = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, true, true);
  }

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1, cvs2;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // we need dg to use correct scaling of basis functions
  Teuchos::ParameterList dglist = plist.sublist("PK operator")
                                       .sublist("flux operator").sublist("schema");
  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dglist, mesh0));

  AnalyticDG04 ana(mesh0, order, true);
  ana.InitialGuess(*dg, p1c, 1.0);

  // create remap object
  MyRemapDG remap(mesh0, mesh1, plist);
  DeformMesh(mesh1, deform, 1.0);
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
  CompositeVector p1aux(*p1);
  Explicit_TI::RK<CompositeVector> rk(remap, rk_method, p1aux);

  remap.NonConservativeToConservative(0.0, *p1, p1aux);

  double t(0.0), tend(1.0);
  while(t < tend - dt/2) {
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

  CompositeVector q2(p2), perr(p2);
  Epetra_MultiVector& q2c = *q2.ViewComponent("cell");
  Epetra_MultiVector& pec = *perr.ViewComponent("cell");
  q2c = p2c;

  double pnorm, l2_err, inf_err, l20_err, inf0_err;
  ana.ComputeCellErrorRemap(*dg, p2c, tend, 0, mesh1,
                            pnorm, l2_err, inf_err, l20_err, inf0_err, &pec);

  CHECK(((dim == 2) ? l2_err : l20_err) < 0.12 / (order + 1));

  if (MyPID == 0) {
    printf("nx=%3d (orig) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n", 
        nx, l20_err, l2_err, inf0_err, inf_err);
  }

  // optional projection on the space of polynomials 
  for (int c = 0; c < ncells_owned; ++c) {
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

  ana.ComputeCellErrorRemap(*dg, q2c, tend, 1, mesh1,
                            pnorm, l2_err, inf_err, l20_err, inf0_err);

  if (MyPID == 0) {
    printf("nx=%3d (proj) L2=%12.8g(mean) %12.8g  Inf=%12.8g %12.8g\n", 
        nx, l20_err, l2_err, inf0_err, inf_err);
  }

  // concervation errors: mass and volume (CGL)
  auto& det = remap.det();
  double area(0.0), area1(0.0), mass1(0.0), gcl_err(0.0), gcl_inf(0.0);
  WhetStone::NumericalIntegration<AmanziMesh::MeshLight> numi(mesh0);

  for (int c = 0; c < ncells_owned; ++c) {
    double vol1 = numi.IntegratePolynomialCell(c, det[c].Value(1.0));
    double vol2 = mesh1->cell_volume(c);

    area += vol1;
    area1 += vol2;

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
  ana.GlobalOp("sum", &area1, 1);
  ana.GlobalOp("sum", &mass1, 1);
  ana.GlobalOp("sum", &gcl_err, 1);
  ana.GlobalOp("max", &gcl_inf, 1);

  if (MyPID == 0) {
    printf("Conservation: dMass=%10.4g  dVol=%10.6g  dVolLinear=%10.6g\n",
           mass1 - mass0, 1.0 - area, 1.0 - area1);
    printf("GCL: L1=%12.8g  Inf=%12.8g\n", gcl_err, gcl_inf);
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh1, true, false);

  io.InitializeCycle(t, 1, "");
  io.WriteVector(*p2c(0), "solution", AmanziMesh::CELL);
  io.WriteVector(*q2c(0), "solution-prj", AmanziMesh::CELL);
  io.FinalizeCycle();
}

TEST(REMAP_DUAL_2D) {
  std::string maps = "VEM";
  double dT(0.1);
  RemapTestsDualRK("FEM", "", 10,10,0, dT);
  RemapTestsDualRK(maps, "test/median15x16.exo", 16,1,0, dT/2);
}

TEST(REMAP_DUAL_3D) {
  std::string maps = "VEM";
  double dT(0.1);
  int deform = 1;
  RemapTestsDualRK(maps, "", 4,4,4, dT, deform);
}

TEST(REMAP_DUAL_DEV) {
  /*
  double dT(0.025);
  int deform = 5;
  RemapTestsDualRK(maps, "test/prism10.exo", 10,1,1, dT,   deform);
  RemapTestsDualRK(maps, "test/prism20.exo", 20,1,1, dT/2, deform);
  RemapTestsDualRK(maps, "test/prism40.exo", 40,1,1, dT/4, deform);
  */

  /*
  double dT(0.025);
  int deform = 5;
  RemapTestsDualRK(maps, "test/hexes4.exo",   4,1,1, dT,   deform);
  RemapTestsDualRK(maps, "test/hexes8.exo",   8,1,1, dT/2, deform);
  RemapTestsDualRK(maps, "test/hexes16.exo", 16,1,1, dT/4, deform);
  RemapTestsDualRK(maps, "test/hexes32.exo", 32,1,1, dT/8, deform);
  */

  /*
  double dT(0.1);
  int deform = 1;
  RemapTestsDualRK(maps, "", 10,10,10, dT/2, deform);
  RemapTestsDualRK(maps, "", 20,20,20, dT/4, deform);
  RemapTestsDualRK(maps, "", 40,40,40, dT/8, deform);
  */

  /*
  double dT(0.01);
  int deform = 1;
  RemapTestsDualRK(maps, "",  16, 16,0, dT,    deform);
  RemapTestsDualRK(maps, "",  32, 32,0, dT/2,  deform);
  RemapTestsDualRK(maps, "",  64, 64,0, dT/4,  deform);
  RemapTestsDualRK(maps, "", 128,128,0, dT/8,  deform);
  RemapTestsDualRK(maps, "", 256,256,0, dT/16, deform);
  */

  /*
  double dT(0.01);
  int deform = 5;
  RemapTestsDualRK(maps, "test/median15x16.exo",    16,0,0, dT,   deform);
  RemapTestsDualRK(maps, "test/median32x33.exo",    32,0,0, dT/2, deform);
  RemapTestsDualRK(maps, "test/median63x64.exo",    64,0,0, dT/4, deform);
  RemapTestsDualRK(maps, "test/median127x128.exo", 128,0,0, dT/8, deform);
  RemapTestsDualRK(maps, "test/median255x256.exo", 256,0,0, dT/16,deform);
  */

  /*
  double dT(0.05);
  int deform = 5;
  RemapTestsDualRK(maps, "test/mesh_poly20x20.exo",    20,0,0, dT,   deform);
  RemapTestsDualRK(maps, "test/mesh_poly40x40.exo",    40,0,0, dT/2, deform);
  RemapTestsDualRK(maps, "test/mesh_poly80x80.exo",    80,0,0, dT/4, deform);
  RemapTestsDualRK(maps, "test/mesh_poly160x160.exo", 160,0,0, dT/8, deform);
  */

  /*
  double dT(0.05);
  int deform = 5;
  RemapTestsDualRK(maps, "test/random10.exo", 10,0,0, dT,   deform);
  RemapTestsDualRK(maps, "test/random20.exo", 20,0,0, dT/2, deform);
  RemapTestsDualRK(maps, "test/random40.exo", 40,0,0, dT/4, deform);
  */

  /*
  double dT(0.025);
  RemapTestsDualRK(maps, "test/triangular8.exo",  0,0,0, dT);
  RemapTestsDualRK(maps, "test/triangular16.exo", 0,0,0, dT/2);
  RemapTestsDualRK(maps, "test/triangular32.exo", 0,0,0, dT/4);
  RemapTestsDualRK(maps, "test/triangular64.exo", 0,0,0, dT/8);
  RemapTestsDualRK(maps, "test/triangular128.exo",0,0,0, dT/16);
  */
}

