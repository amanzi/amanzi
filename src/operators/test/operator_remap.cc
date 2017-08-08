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
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshMaps_FEM.hh"
#include "MeshMaps_VEM.hh"
#include "Tensor.hh"

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
  // Teuchos::RCP<const Mesh> mesh0 = meshfactory("test/median32x33.exo", Teuchos::null);

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_wghost = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_owned = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  // create second and auxiliary mesh
  Teuchos::RCP<Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
  // Teuchos::RCP<Mesh> mesh1 = meshfactory("test/median32x33.exo", Teuchos::null);

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
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

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

  // create advection operator
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

  Teuchos::RCP<std::vector<WhetStone::Polynomial> > jac0 = 
     Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_owned));
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > jac1 = 
     Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_owned));

  op_reac0->Setup(jac0);
  op_reac1->Setup(jac1);

  double t(0.0), tend(1.0);
  WhetStone::Polynomial det0, det1;
  WhetStone::MeshMaps_FEM maps(mesh0, mesh1);

  while(t < tend - dt/2) {
    // calculate determinant of Jacobian at time t
    for (int c = 0; c < ncells_owned; ++c) {
      maps.JacobianDet(c, t, (*jac0)[c]);
      maps.JacobianDet(c, t + dt, (*jac1)[c]);
    }

    // rotate velocities and calculate normal component
    Entity_ID_List faces;
    std::vector<int> dirs;
    WhetStone::Polynomial poly(2, 1);
    std::vector<WhetStone::Polynomial> uv(2, poly);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh0)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);

    Teuchos::RCP<CompositeVector> vel = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs));
    Epetra_MultiVector& vel_f = *vel->ViewComponent("face");

    for (int c = 0; c < ncells_owned; ++c) {
      mesh0->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      double sum0(0.0), sum1(0.0);
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];

        // calculate j J^{-t} N dA
        maps.VelocityFace(c, f, uv);

        WhetStone::Tensor J(2, 2); 
        maps.JacobianFaceValue(c, f, uv, xref, J);
        J *= t;
        J += 1.0 - t;
        WhetStone::Tensor C = J.Cofactors();
        AmanziGeometry::Point cn = C * mesh0->face_normal(f); 

        // calculate velocity
        const AmanziGeometry::Point& xf = mesh0->face_centroid(f);
        xv = AmanziGeometry::Point(uv[0].Value(xf), uv[1].Value(xf));
        vel_f[0][f] = xv * cn;

        // test
        /*
        const AmanziGeometry::Point& xf1 = mesh1->face_centroid(f);
        sum0 += (xf + t * (xf1 - xf)) * cn * dirs[n];
        maps.JacobianFaceValue(c, f, uv, xref, J);
        J *= t + dt;
        J += 1.0 - (t + dt);
        C = J.Cofactors();
        cn = C * mesh0->face_normal(f); 
        sum1 += (xf + (t + dt) * (xf1 - xf)) * cn * dirs[n];
        */
      }
      // (*jac0)[0][c] = sum0 / 2 / mesh0->cell_volume(c);
      // (*jac1)[0][c] = sum1 / 2 / mesh0->cell_volume(c);
    }

    // populate operators
    op->UpdateMatrices(vel.ptr());
    op_reac0->UpdateMatrices(p1.ptr());

    // predictor step
    CompositeVector& rhs = *global_reac0->rhs();
    global_reac0->Apply(*p1, rhs);

    CompositeVector g(cvs1);
    global_op->Apply(*p1, g);
    g.Update(1.0, rhs, dt);

    op_reac1->UpdateMatrices(p1.ptr());
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

    *p1->ViewComponent("cell") = *p2.ViewComponent("cell");
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
  RemapTests2DExplicit(0, "DG order 0", 30, 30, 0.1 / 3);
}

// TEST(REMAP_DG1_EXPLICIT) {
//   RemapTests2DExplicit(1, "DG order 1", 20, 20, 0.1 / 2);
// }

