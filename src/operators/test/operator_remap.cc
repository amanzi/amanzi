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

#include "UnitTest++.h"

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "LinearOperatorPCG.hh"

#include "Accumulation.hh"
#include "AdvectionRiemann.hh"
#include "OperatorDefs.hh"
#include "Reaction.hh"
#include "RemapUtils.hh"


/* *****************************************************************
* Remap of polynomilas in two dimensions
***************************************************************** */
void RemapTests2D(int order, std::string disc_name) {
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

  int nx(10), ny(10);
  // Teuchos::RCP<const Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
  Teuchos::RCP<const Mesh> mesh1 = meshfactory("test/random10.exo");

  int ncells_owned = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_wghost = mesh1->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  // create deformed mesh
  // Teuchos::RCP<Mesh> mesh2 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
  Teuchos::RCP<Mesh> mesh2 = meshfactory("test/random10.exo");

  int nnodes_owned = mesh2->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  AmanziGeometry::Point xv(2);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh2->node_get_coordinates(v, &xv);
    if (!(fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
          fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6)) {
      xv[0] += 0.5 / nx * ((double)rand() / RAND_MAX - 0.5);
      xv[1] += 0.5 / ny * ((double)rand() / RAND_MAX - 0.5);
      nodeids.push_back(v);
      new_positions.push_back(xv);
    }
  }
    
  mesh2->deform(nodeids, new_positions, false, &final_positions);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p1(cvs1);
  Epetra_MultiVector& p1c = *p1.ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh1->cell_centroid(c);
    p1c[0][c] = xc[0] * xc[0] + 2 * xc[1];
    if (nk > 1) {
      p1c[1][c] = 2.0 * xc[0];
      p1c[2][c] = 2.0;
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
  plist.set<std::string>("discretization", disc_name);

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

  // calculate mesh velocity on faces and in cells
  Teuchos::RCP<CompositeVector> velc, velf;

  RemapVelocityFaces(order, mesh1, mesh2, velf);
  RemapVelocityCells(order, mesh1, mesh2, velc);

  // populate operators
  op->UpdateMatrices(*velf);
  op_adv->UpdateMatrices(*velc);
  op_reac->UpdateMatrices(p1);

  // predictor step
  CompositeVector& rhs = *global_reac->rhs();
  global_reac->Apply(p1, rhs);

  CompositeVector g(cvs1);
  global_op->Apply(p1, g);
  g.Update(1.0, rhs, 1.0);

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

  // calculate error
  double pl2_err(0.0), pinf_err(0.0), area(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh2->cell_centroid(c);
    double tmp = (xc[0] * xc[0] + 2 * xc[1]) - p2c[0][c];

    pinf_err = std::max(pinf_err, fabs(tmp));
    pl2_err += tmp * tmp * mesh2->cell_volume(c);
    area += mesh2->cell_volume(c);
  }
  pl2_err = std::pow(pl2_err, 0.5);

  if (MyPID == 0) {
    printf("L2(p)=%12.8g  Inf(p)=%12.8g  Err(area)=%12.8g\n", 
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


TEST(REMAP_2D) {
  RemapTests2D(1, "DG order 1");
}

