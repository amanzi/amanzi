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
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "OutputXDMF.hh"

// Amanzi::Operators
#include "MeshDeformation.hh"
#include "RemapDG_Tests.hh"

#include "AnalyticDG04.hh"

namespace Amanzi {

class MyRemapDG : public RemapDG_Tests<AnalyticDG04> {
 public:
  MyRemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
            const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
            Teuchos::ParameterList& plist)
    : RemapDG_Tests<AnalyticDG04>(mesh0, mesh1, plist) {};
  ~MyRemapDG() {};

  // access 
  const std::vector<WhetStone::SpaceTimePolynomial> det() const { return *det_; }
};


/* *****************************************************************
* Initialization of the consistent jacobian determinant
***************************************************************** */
template<class AnalyticDG>
void MyRemapDG<AnalyticDG>::InitializeConsistentJacobianDeterminant()
{
  // constant part of determinant
  op_adv_->Setup(velc_, false);
  op_adv_->UpdateMatrices(0.0);

  op_reac_->Setup(det_, false);
  op_reac_->UpdateMatrices(0.0);

  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

  op_flux_->Setup(velf_.ptr(), false);
  op_flux_->UpdateMatrices(0.0);
  op_flux_->ApplyBCs(true, true, true);

  CompositeVector& tmp = *op_reac_->global_operator()->rhs();
  CompositeVector one(tmp), u0(tmp), u1(tmp);
  Epetra_MultiVector& one_c = *one.ViewComponent("cell", true);

  one.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_wghost_; ++c) one_c[0][c] = 1.0;

  op_flux_->global_operator()->Apply(one, tmp);
  op_reac_->global_operator()->Apply(tmp, u0);

  // linear part of determinant
  double dt(0.01);
  op_adv_->UpdateMatrices(dt);

  op_flux_->UpdateMatrices(dt);
  op_flux_->ApplyBCs(true, true, true);

  op_flux_->global_operator()->Apply(one, tmp);
  op_reac_->global_operator()->Apply(tmp, u1);
  u1.Update(-1.0/dt, u0, 1.0/dt);

  // save as polynomials
  int nk = one_c.NumVectors();
  Amanzi::WhetStone::DenseVector data(nk);
  Epetra_MultiVector& u0c = *u0.ViewComponent("cell", true);
  Epetra_MultiVector& u1c = *u1.ViewComponent("cell", true);

  det0_.resize(ncells_owned_);
  det1_.resize(ncells_owned_);

  for (int c = 0; c < ncells_owned_; ++c) {
    const auto& basis = dg_->cell_basis(c);

    for (int i = 0; i < nk; ++i) data(i) = u0c[i][c];
    det0_[c] = basis.CalculatePolynomial(mesh0_, c, order_, data);

    for (int i = 0; i < nk; ++i) data(i) = u1c[i][c];
    det1_[c] = basis.CalculatePolynomial(mesh0_, c, order_, data);
  }
}

}  // namespace Amanzi


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapGCL(const Amanzi::Explicit_TI::method_t& rk_method,
              std::string map_name, std::string file_name,
              int nx, int ny, int nz, double dt,
              int deform = 1) {
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
  plist.sublist("maps").set<std::string>("map name", map_name);

  // print simulation header
  const auto& map_list = plist.sublist("maps");
  int vel_order = map_list.get<int>("method order");

  if (MyPID == 0) {
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
    std::string map_name = map_list.get<std::string>("map name");
      
    std::cout << "\nTest: " << dim << "D remap:"
              << " mesh=" << ((file_name == "") ? "structured" : file_name)
              << " deform=" << deform << std::endl;

    std::cout << "      discretization: order=" << order 
              << ", map=" << map_name << std::endl;

    std::cout << "      map details: order=" << vel_order 
              << ", projector=" << vel_projector 
              << ", method=\"" << vel_method << "\"" << std::endl;
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
  if (MyPID == 0) std::cout << "Deforming mesh...\n";
  DeformMesh(mesh1, deform, 1.0);
  remap.InitializeOperators(dg);
  if (MyPID == 0) std::cout << "Computing static data on mesh scheleton...\n";
  remap.StaticEdgeFaceVelocities();
  remap.StaticFaceCoVelocity();
  if (MyPID == 0) std::cout << "Computing static data in mesh cells...\n";
  remap.StaticCellVelocity();
  remap.StaticCellCoVelocity();
  if (MyPID == 0) std::cout << "Done.\n";

  // initial mass
  double mass0(0.0);
  WhetStone::DenseVector data(nk);
  WhetStone::NumericalIntegration<AmanziMesh::Mesh> numi(mesh0);

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

  CompositeVectorSpace cvs3;
  cvs3.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

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

  // concervation errors: mass and volume (CGL)
  auto& det = remap.det();
  double area(0.0), area1(0.0), mass1(0.0), gcl_err(0.0), gcl_inf(0.0);

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

    int quad_order = det[c].order() + poly.order();

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

  io.InitializeCycle(t, 1);
  io.WriteVector(*p2c(0), "solution", AmanziMesh::CELL);
  io.WriteVector(*q2c(0), "solution-prj", AmanziMesh::CELL);
  io.FinalizeCycle();
}

TEST(REMAP_GEOMETRIC_CONSERVATION_LAW) {
  double dT(0.1);
  auto rk_method = Amanzi::Explicit_TI::heun_euler;
  int deform = 1;
  RemapGCL(rk_method, "VEM", "test/median15x16.exo", 16,1,0, dT/2);
}

