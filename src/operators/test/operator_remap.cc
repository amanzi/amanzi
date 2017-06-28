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
    p1c[0][c] = xc[0];  // + 2 * xc[1] * xc[1];
    if (nk > 1) {
      p1c[1][c] = 1.0;
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

  // create secondary advection operator (cell-based)
  plist.sublist("schema domain")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op_adv = Teuchos::rcp(new AdvectionRiemann(plist, global_op));

  // create accumulation operator
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<Reaction> op_reac = Teuchos::rcp(new Reaction(plist, mesh1));
  auto global_reac = op_reac->global_operator();

  Teuchos::RCP<Epetra_MultiVector> jac = Teuchos::rcp(new Epetra_MultiVector(mesh1->cell_map(true), 1));
  op_reac->Setup(jac);

  double t(0.0), tend(1.0);
  while(t < tend) {
    // calculate determinat of jacobian
    for (int c = 0; c < ncells_owned; ++c) {
      double v0 = mesh1->cell_volume(c);
      double v1 = mesh2->cell_volume(c);
      (*jac)[0][c] = 1.0 + t * (v1 - v0) / v0;
    }

    // rotate velocities
    WhetStone::Tensor R(2, 2);
    CompositeVector velc_t(*velc), velf_t(*velf);

    Epetra_MultiVector& vel = *velf_t.ViewComponent("face");
double sum(0.0);
    for (int f = 0; f < nfaces_wghost; ++f) {
      R(0, 0) = 1.0 + t * vel[5][f];
      R(0, 1) = -t * vel[4][f];
      R(1, 0) = -t * vel[3][f];
      R(1, 1) = 1.0 + t * vel[2][f];

      for (int k = 0; k < nk; ++k) {
        xv[0] = vel[2 * k][f];
        xv[1] = vel[2 *k + 1][f];
if (f < 4) {
sum += xv * mesh1->face_normal(f);
AmanziGeometry::Point xv2(2), xv3(2);
xv2[0] = vel[2][f];
xv2[1] = vel[3][f];
xv3[0] = vel[4][f];
xv3[1] = vel[5][f];
std::cout << f << " " << xv << " " << xv2 << " " << xv3 << " sum=" << sum << " nf=" << mesh1->face_normal(f) << std::endl;
}

        xv = R * xv;
      
        vel[2 * k][f] = xv[0];
        vel[2 * k + 1][f] = xv[1];
      }
    }
    
    Epetra_MultiVector& vel2 = *velc_t.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; ++c) {
      R(0, 0) = 1.0 + t * vel2[5][c];
      R(0, 1) = -t * vel2[4][c];
      R(1, 0) = -t * vel2[3][c];
      R(1, 1) = 1.0 + t * vel2[2][c];

      for (int k = 0; k < nk; ++k) {
        xv[0] = vel2[2 * k][c];
        xv[1] = vel2[2 *k + 1][c];

        xv = R * xv;
      
        vel2[2 * k][c] = xv[0];
        vel2[2 * k + 1][c] = xv[1];
      }
    }
    
    // populate operators
    op->UpdateMatrices(velf_t);
    op_adv->UpdateMatrices(velc_t);
    op_reac->UpdateMatrices(p1);

    // predictor step
    CompositeVector& rhs = *global_reac->rhs();
    global_reac->Apply(p1, rhs);

    CompositeVector g(cvs1);
    global_op->Apply(p1, g);
    g.Update(1.0, rhs, dt);

    global_reac->SymbolicAssembleMatrix();
    global_reac->AssembleMatrix();

    plist.set<std::string>("preconditioner type", "diagonal");
    global_reac->InitPreconditioner(plist);

    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        pcg(global_reac, global_reac);

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

    double tmp = (xc[0] + 0 * xc[1] * xc[1]) - p2c[0][c];
    pinf_err = std::max(pinf_err, fabs(tmp));
    pl2_err += tmp * tmp * area_c;
    // std::cout << c << " " << p2c[0][c] << " err=" << tmp << std::endl;

    area += area_c;
  }
  pl2_err = std::pow(pl2_err, 0.5);

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


/* *****************************************************************
* Remap of polynomilas in three dimensions. Implicit scheme.
***************************************************************** */
void RemapTests2DImplicit(int order, std::string disc_name,
                          int nx, int ny, double dt) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: remap of functions in 2D: implicit" << std::endl;

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
    p1c[0][c] = xc[0];  // + 2 * xc[1] * xc[1];
    if (nk > 1) {
      p1c[1][c] = 1.0;
      p1c[2][c] = 0.0;  // 4.0 * xc[1];
    }
  }

  // allocate memory
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh2)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell", true);

  // create accumulation operator
  Teuchos::ParameterList plist;
  plist.set<std::string>("discretization", disc_name);

  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<Reaction> op_reac = Teuchos::rcp(new Reaction(plist, mesh1));
  auto global_reac = op_reac->global_operator();

  // create primary advection operator
  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op = Teuchos::rcp(new AdvectionRiemann(plist, mesh1));
  auto global_op = op->global_operator();

  // add secondary advection operator (cell-based)
  plist.sublist("schema domain")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op_adv = Teuchos::rcp(new AdvectionRiemann(plist, global_op));

  // add accumulation operator
  global_op->OpPushBack(op_reac->local_operator());

  // START TIME INTEGRATION
  Teuchos::RCP<Epetra_MultiVector> jac = Teuchos::rcp(new Epetra_MultiVector(mesh1->cell_map(true), 1));
  op_reac->Setup(jac);

  double t(0.0), tend(1.0);
  while(t < tend) {
    t += dt;

    // calculate determinant of Jacobian matrix
    for (int c = 0; c < ncells_owned; ++c) {
      double v0 = mesh1->cell_volume(c);
      double v1 = mesh2->cell_volume(c);
      (*jac)[0][c] = 1.0 + t * (v1 - v0) / v0;
    }

    // rotate velocities
    WhetStone::Tensor R(2, 2);
    CompositeVector velc_t(*velc), velf_t(*velf);

    Epetra_MultiVector& vel = *velf_t.ViewComponent("face");
    for (int f = 0; f < nfaces_wghost; ++f) {
      R(0, 0) = 1.0 + t * vel[5][f];
      R(0, 1) = -t * vel[4][f];
      R(1, 0) = -t * vel[3][f];
      R(1, 1) = 1.0 + t * vel[2][f];

      for (int k = 0; k < nk; ++k) {
        xv[0] = vel[2 * k][f];
        xv[1] = vel[2 *k + 1][f];

        xv = R * xv;
      
        vel[2 * k][f] = xv[0];
        vel[2 * k + 1][f] = xv[1];
      }
    }
    
    Epetra_MultiVector& vel2 = *velc_t.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; ++c) {
      R(0, 0) = 1.0 + t * vel2[5][c];
      R(0, 1) = -t * vel2[4][c];
      R(1, 0) = -t * vel2[3][c];
      R(1, 1) = 1.0 + t * vel2[2][c];

      for (int k = 0; k < nk; ++k) {
        xv[0] = vel2[2 * k][c];
        xv[1] = vel2[2 *k + 1][c];

        xv = R * xv;
      
        vel2[2 * k][c] = xv[0];
        vel2[2 * k + 1][c] = xv[1];
      }
    }
    
    // populate operators
    op->UpdateMatrices(velf_t);
    op_adv->UpdateMatrices(velc_t);
    op_reac->UpdateMatrices(p1);

    op->local_operator()->Rescale(-dt);
    op_adv->local_operator()->Rescale(-dt);

    // predictor step
    CompositeVector& rhs = *global_reac->rhs();
    global_reac->Apply(p1, rhs);

    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();

    plist.set<std::string>("preconditioner type", "diagonal");
    global_op->InitPreconditioner(plist);

    AmanziSolvers::LinearOperatorGMRES<Operator, CompositeVector, CompositeVectorSpace>
        gmres(global_op, global_op);

    plist.set<double>("error tolerance", 1e-10)
         .set<int>("maximum number of iterations", 300);

    gmres.Init(plist);
    gmres.ApplyInverse(rhs, p2);

    *p1.ViewComponent("cell") = *p2.ViewComponent("cell");
  }

  // calculate error
  double pl2_err(0.0), pinf_err(0.0), area(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh2->cell_centroid(c);
    double area_c = mesh2->cell_volume(c);

    double tmp = (xc[0] + 0 * xc[1] * xc[1]) - p2c[0][c];
    pinf_err = std::max(pinf_err, fabs(tmp));
    pl2_err += tmp * tmp * area_c;

    area += area_c;
  }
  pl2_err = std::pow(pl2_err, 0.5);

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
  RemapTests2DExplicit(0, "DG order 0", 40, 40, 0.1);
}

// TEST(REMAP_2D_IMPLICIT) {
//   RemapTests2DImplicit(0, "DG order 0", 40, 40, 0.1);
// }

