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
#include "MeshFactory.hh"

#include "Accumulation.hh"
#include "AdvectionRiemann.hh"
#include "OperatorDefs.hh"


/* *****************************************************************
* Remap of polynomilas in two dimensions
***************************************************************** */
TEST(REMAP_CONSTANT_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: remap of constant functions in 2D." << std::endl;

  // create initial mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  int nx(32), ny(32);
  Teuchos::RCP<const Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);

  int ncells_owned = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_wghost = mesh1->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  // create deformed mesh
  Teuchos::RCP<Mesh> mesh2 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);

  int nnodes_owned = mesh2->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  AmanziGeometry::Point xv(2);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh2->node_get_coordinates(v, &xv);
    if (!(fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
          fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6)) {
      xv[0] += 0.02 / 4;
      xv[1] += 0.02 / 4;
      nodeids.push_back(v);
      new_positions.push_back(xv);
    }
  }
    
  mesh2->deform(nodeids, new_positions, false, &final_positions);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);
  CompositeVector p1(cvs1);
  Epetra_MultiVector& p1c = *p1.ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh1->cell_centroid(c);
    p1c[0][c] = xc[0] + 2 * xc[1];
  }

  // remap cell-based field
  // -- allocate memory
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh2)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell", true);

  // -- create primary advection operator
  Teuchos::ParameterList plist;
  plist.set<std::string>("discretization", "DG order 0: face");

  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({1}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op = Teuchos::rcp(new AdvectionRiemann(plist, mesh1));
  auto global_op = op->global_operator();

  // -- create secondary advection operator (cell-based)
  plist.set<std::string>("discretization", "DG order 0: cell");

  plist.sublist("schema domain")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({1}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op_adv = Teuchos::rcp(new AdvectionRiemann(plist, global_op));

  // -- create accumulation operator
  Teuchos::RCP<Accumulation> op_acc = Teuchos::rcp(new Accumulation(AmanziMesh::CELL, mesh1));
  auto global_acc = op_acc->global_operator();

  // -- calculate flux on mesh faces
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh1)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);
  CompositeVector flux(cvs);
  Epetra_MultiVector& flux_f = *flux.ViewComponent("face", true);

  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf1 = mesh1->face_centroid(f);
    const AmanziGeometry::Point& xf2 = mesh2->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh1->face_normal(f);
    double area = mesh1->face_area(f);

    flux_f[0][f] = (xf2 - xf1) * normal / area;
  }

  // -- populate operators
  op->UpdateMatrices(flux);
  op_acc->AddAccumulationDelta(p1, 1.0, "cell");

  // -- invert the mass matrix
  CompositeVector g(cvs1);
  global_op->Apply(p1, g);

  CompositeVector& rhs = *global_acc->rhs();
  g.Update(1.0, rhs, -1.0);

  global_acc->SymbolicAssembleMatrix();
  global_acc->AssembleMatrix();

  plist.set<std::string>("preconditioner type", "diagonal");
  global_acc->InitPreconditioner(plist);
  global_acc->ApplyInverse(g, p2);

  // calculate error
  double pl2_err(0.0), pinf_err(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh2->cell_centroid(c);
    double tmp = (xc[0] + 2 * xc[1]) - p2c[0][c];

    pinf_err = std::max(pinf_err, fabs(tmp));
    pl2_err += tmp * tmp * mesh2->cell_volume(c);
  }
  pl2_err = std::pow(pl2_err, 0.5);

  if (MyPID == 0) {
    printf("L2(p)=%12.8g  Inf(p)=%12.8g\n", pl2_err, pinf_err);
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



