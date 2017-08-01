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

#include "DG_Modal.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

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

  Teuchos::RCP<const Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);

  int ncells_owned = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_wghost = mesh1->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_owned = mesh1->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  // create second and auxiliary mesh
  Teuchos::RCP<Mesh> mesh2 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);

  // deform the second mesh
  AmanziGeometry::Point xv(2);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh2->node_get_coordinates(v, &xv);

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
  mesh2->deform(nodeids, new_positions, false, &final_positions);

  // calculate mesh velocity on faces and in cells
  Teuchos::RCP<CompositeVector> velc, velf;
  RemapVelocityFaces(1, mesh1, mesh2, velf);
  RemapVelocityCells(1, mesh1, mesh2, velc);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p1(cvs1);
  Epetra_MultiVector& p1c = *p1.ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh1->cell_centroid(c);
    p1c[0][c] = xc[0] + 2 * xc[1] * xc[1];
    if (nk > 1) {
      p1c[1][c] = 0.0;
      p1c[2][c] = 0.0;  // 4.0 * xc[1];
    }
  }

  // allocate memory
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh2)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
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

  Teuchos::RCP<AdvectionRiemann> op = Teuchos::rcp(new AdvectionRiemann(plist, mesh1));
  auto global_op = op->global_operator();

  // create accumulation operator
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<Reaction> op_reac0 = Teuchos::rcp(new Reaction(plist, mesh1));
  Teuchos::RCP<Reaction> op_reac1 = Teuchos::rcp(new Reaction(plist, mesh1));
  auto global_reac0 = op_reac0->global_operator();
  auto global_reac1 = op_reac1->global_operator();

  Teuchos::RCP<Epetra_MultiVector> jac0 = Teuchos::rcp(new Epetra_MultiVector(mesh1->cell_map(true), 1));
  Teuchos::RCP<Epetra_MultiVector> jac1 = Teuchos::rcp(new Epetra_MultiVector(mesh1->cell_map(true), 1));
  op_reac0->Setup(jac0);
  op_reac1->Setup(jac1);

  double t(0.0), tend(1.0);
  WhetStone::DG_Modal dg1(mesh1), dg2(mesh1, mesh2);

  while(t < tend - dt/2) {
    // calculate determinant of Jacobian
    for (int c = 0; c < ncells_owned; ++c) {
      AmanziGeometry::Point xref(0.5, 0.5);

      WhetStone::Tensor J0 = dg2.EvaluateJacobian(c, xref);
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
    std::vector<int> dirs;

    CompositeVector velf_tmp(*velf);
    const Epetra_MultiVector& vel = *velf->ViewComponent("face");
    Epetra_MultiVector& vel_tmp = *velf_tmp.ViewComponent("face");

    for (int c = 0; c < ncells_owned; ++c) {
      mesh1->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];

        // calculate j J^{-t} N dA
        // double qa(WhetStone::q1d_points[1][0]), qb(WhetStone::q1d_points[1][1]);
        double qa(0.5);
        std::vector<AmanziGeometry::Point> xref;
        if (n == 0) { 
          xref.push_back(AmanziGeometry::Point(qa, 0.0));
          // xref.push_back(AmanziGeometry::Point(qb, 0.0));
        } else if (n == 1) {
          xref.push_back(AmanziGeometry::Point(1.0, qa));
          // xref.push_back(AmanziGeometry::Point(1.0, qb));
        } else if (n == 2) { 
          xref.push_back(AmanziGeometry::Point(qa, 1.0));
          // xref.push_back(AmanziGeometry::Point(qb, 1.0));
        } else if (n == 3) {
          xref.push_back(AmanziGeometry::Point(0.0, qa));
          // xref.push_back(AmanziGeometry::Point(0.0, qb));
        }

        vel_tmp[0][f] = 0.0;
        for (int i = 0; i < xref.size(); ++i) {
          WhetStone::Tensor J = dg2.EvaluateJacobian(c, xref[i]);
          J *= t;
          J += 1.0 - t;
          WhetStone::Tensor C = J.Cofactors();
          AmanziGeometry::Point cn = C * mesh1->face_normal(f); 

          // calculate velocity
          xv = dg2.EvaluateMap(c, xref[i]) - dg1.EvaluateMap(c, xref[i]);
          vel_tmp[0][f] += (xv * cn) * WhetStone::q1d_weights[0][i];
        }
      }
    }

    // populate operators
    op->UpdateMatrices(velf_tmp);
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
    const AmanziGeometry::Point& xc = mesh2->cell_centroid(c);
    double area_c = mesh2->cell_volume(c);

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
    GMV::open_data_file(*mesh2, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p2c, 0, "remaped");
    GMV::close_data_file();
  }
}


TEST(REMAP_2D_EXPLICIT) {
  RemapTests2DExplicit(0, "DG order 0", 20, 20, 0.1 / 2);
}

// TEST(REMAP_2D_IMPLICIT) {
//   RemapTests2DImplicit(0, "DG order 0", 40, 40, 0.1);
// }

