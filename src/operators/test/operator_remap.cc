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

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "DG_Modal.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Accumulation.hh"
#include "AdvectionRiemann.hh"
#include "OperatorDefs.hh"
#include "Reaction.hh"
#include "RemapUtils.hh"


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
***************************************************************** */
void RemapTests2DExplicit(int order, std::string disc_name,
                          int nx, int ny, double dt) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: remap of constant functions in 2D." << std::endl;

  // polynomial space
  int nk = (order + 1) * (order + 2) / 2;

  // create initial mesh
  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK}));

  Teuchos::RCP<const Mesh> mesh0 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_wghost = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_owned = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  // create second and auxiliary mesh
  Teuchos::RCP<Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);

  // deform the second mesh
  AmanziGeometry::Point xv(2), xref(2);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    double ds(0.001), ux, uy;
    for (int i = 0; i < 1000; ++i) {
      ux = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
      uy =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);

      xv[0] += ux * ds;
      xv[1] += uy * ds;
    }

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p1(cvs1);
  Epetra_MultiVector& p1c = *p1.ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh0->cell_centroid(c);
    p1c[0][c] = xc[0] + 2 * xc[1] * xc[1];
    if (nk > 1) {
      p1c[1][c] = 1.0;
      p1c[2][c] = 4.0 * xc[1];
    }
  }

  // allocate memory
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell", true);

  // create primary advection operator
  Teuchos::ParameterList plist;
  plist.set<std::string>("discretization", disc_name);

  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op = Teuchos::rcp(new AdvectionRiemann(plist, mesh0));
  auto global_op = op->global_operator();

  // create accumulation operator
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<Reaction> op_reac0 = Teuchos::rcp(new Reaction(plist, mesh0));
  Teuchos::RCP<Reaction> op_reac1 = Teuchos::rcp(new Reaction(plist, mesh0));
  auto global_reac0 = op_reac0->global_operator();
  auto global_reac1 = op_reac1->global_operator();

  Teuchos::RCP<Epetra_MultiVector> jac0 = Teuchos::rcp(new Epetra_MultiVector(mesh0->cell_map(true), 1));
  Teuchos::RCP<Epetra_MultiVector> jac1 = Teuchos::rcp(new Epetra_MultiVector(mesh0->cell_map(true), 1));
  op_reac0->Setup(jac0);
  op_reac1->Setup(jac1);

  double t(0.0), tend(1.0);
  WhetStone::DG_Modal dg(mesh0, mesh1);
  dg.set_method(WhetStone::WHETSTONE_METHOD_VEM);

  while(t < tend - dt/2) {
    // calculate determinant of Jacobian
    for (int c = 0; c < ncells_owned; ++c) {
      xref.set(0.5, 0.5);

      WhetStone::Tensor J0 = dg.FEM_Jacobian(c, xref);
      WhetStone::Tensor J1(J0);
      J0 *= t;
      J0 += 1.0 - t;
      (*jac0)[0][c] = J0.Det();

      J1 *= t + dt;
      J1 += 1.0 - (t + dt);
      (*jac1)[0][c] = J1.Det();
    }

    // rotate velocities and calculate normal component
    Entity_ID_List faces;
    WhetStone::Polynomial poly(2, 1);
    std::vector<WhetStone::Polynomial> uv(2, poly);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh0)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);

    Teuchos::RCP<CompositeVector> vel = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs));
    Epetra_MultiVector& vel_f = *vel->ViewComponent("face");

    for (int c = 0; c < ncells_owned; ++c) {
      mesh0->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];

        // calculate j J^{-t} N dA
        dg.FaceVelocity(c, f, uv);

        WhetStone::Tensor J = dg.FaceJacobian(c, f, uv, xref);
        J *= t;
        J += 1.0 - t;
        WhetStone::Tensor C = J.Cofactors();
        AmanziGeometry::Point cn = C * mesh0->face_normal(f); 

        // calculate velocity
        const AmanziGeometry::Point& xf = mesh0->face_centroid(f);
        xv = AmanziGeometry::Point(uv[0].Value(xf), uv[1].Value(xf));
        vel_f[0][f] = xv * cn;
      }
    }

    // populate operators
    op->UpdateMatrices(*vel);
    op_reac0->UpdateMatrices(p1);

    // predictor step
    CompositeVector& rhs = *global_reac0->rhs();
    global_reac0->Apply(p1, rhs);

    CompositeVector g(cvs1);
    global_op->Apply(p1, g);
    g.Update(1.0, rhs, dt);

    op_reac1->UpdateMatrices(p1);
    global_reac1->SymbolicAssembleMatrix();
    global_reac1->AssembleMatrix();

    plist.set<std::string>("preconditioner type", "diagonal");
    global_reac1->InitPreconditioner(plist);

    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        pcg(global_reac1, global_reac1);

    pcg.Init(plist);
    pcg.ApplyInverse(g, p2);

    // corrector step
    /*
    p2.Update(0.5, p1, 0.5);
    global_op->Apply(p2, g);
    g.Update(1.0, rhs, 1.0);

    pcg.ApplyInverse(g, p2);
    */

    *p1.ViewComponent("cell") = *p2.ViewComponent("cell");
    t += dt;
  }

  // calculate error
  double pl2_err(0.0), pinf_err(0.0), area(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh1->cell_centroid(c);
    double area_c = mesh1->cell_volume(c);

    double tmp = (xc[0] + 2 * xc[1] * xc[1]) - p2c[0][c];
    pinf_err = std::max(pinf_err, fabs(tmp));
    pl2_err += tmp * tmp * area_c;

    area += area_c;
  }
  pl2_err = std::pow(pl2_err, 0.5);
  CHECK(pl2_err < 0.05);

  if (MyPID == 0) {
    printf("L2(p0)=%12.8g  Inf(p0)=%12.8g  Err(area)=%12.8g\n", 
        pl2_err, pinf_err, 1.0 - area);
  }

  // visualization
  if (MyPID == 0) {
    const Epetra_MultiVector& p2c = *p2.ViewComponent("cell");
    GMV::open_data_file(*mesh1, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p1c, 0, "remaped");
    GMV::close_data_file();
  }
}


TEST(REMAP_DG0_EXPLICIT) {
  RemapTests2DExplicit(0, "DG order 0", 20, 20, 0.1 / 2);
}

// TEST(REMAP_DG1_EXPLICIT) {
//   RemapTests2DExplicit(1, "DG order 1", 20, 20, 0.1 / 2);
// }

