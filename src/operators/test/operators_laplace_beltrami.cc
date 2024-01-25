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
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"
#include "Verification.hh"


void LaplaceBeltramiFlat(std::vector<std::string> surfaces, std::string diff_op)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) {
    std::cout << "\nTest: Laplace Beltrami solver: ";
    for (int i = 0; i < surfaces.size(); ++i)
      std::cout << "\"" << surfaces[i] << "\", ";
    std::cout << diff_op << std::endl;
  }

  // read parameter list
  std::string xmlFileName = "test/operator_laplace_beltrami.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh
  ParameterList region_list = plist.sublist("Regions Flat");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  // extract a manifold mesh
  RCP<const Mesh> surfmesh = meshfactory.create(mesh, surfaces, AmanziMesh::FACE);

  int ncells_owned = surfmesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  int ncells_wghost = surfmesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
  int nfaces_wghost = surfmesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);

  std::cout << "pid=" << MyPID << " cells: " << ncells_owned << " " << ncells_wghost << std::endl;

  // verify one-to-one map (2D-cell -> 3D-face)
  for (int c = 0; c < ncells_wghost; ++c) {
    int g = surfmesh->getEntityParent(AmanziMesh::CELL, c);
    double diff = AmanziGeometry::norm(surfmesh->getCellCentroid(c) - mesh->getFaceCentroid(g));
    CHECK_CLOSE(0.0, diff, 1e-14); 
  }

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    const Point& xc = mesh->getCellCentroid(c);
    Kc(0, 0) = 1.0 + xc[0] * xc[0];
    K->push_back(Kc);
  }

  // create boundary data (no mixed bc)
  Entity_ID_List cells;
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    surfmesh->face_get_cells(f, Parallel_kind::ALL, &cells);
    if (cells.size() == 2) continue;

    const Point& xf = surfmesh->getFaceCentroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = AmanziGeometry::L22(xf);
    }
  }

  // create diffusion operator 
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist(diff_op);
  auto op = Teuchos::rcp(new PDE_DiffusionMFD(olist, surfmesh));
  op->Init();
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  op->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD();
  ver.CheckPreconditionerSPD();

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers").sublist("PCG").sublist("pcg parameters");
  auto solver = Teuchos::rcp(new AmanziSolvers::LinearOperatorPCG<
      Operator, CompositeVector, CompositeVectorSpace>(global_op, global_op));
  solver->Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);
  if (diff_op == "diffusion operator") 
      ver.CheckResidual(solution, 1.0e-12);

  int num_itrs = solver->num_itrs();
  CHECK(num_itrs < 10);
 
  // check bounds of cell-based solution
  const Epetra_MultiVector& p = *solution.viewComponent("cell");
  for (int c = 0; c < p.MyLength(); ++c) {
    CHECK(p[0][c] > 0.0 && p[0][c] < 3.0);
  } 

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver->residual() << " itr=" << num_itrs
              << " code=" << solver->returned_code() << std::endl;

    // visualization
    GMV::open_data_file(*surfmesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}

TEST(LAPLACE_BELTRAMI_FLAT) {
  // boundary surface
  std::vector<std::string> surfaces(1);
  surfaces[0] = "Top surface";
  LaplaceBeltramiFlat(surfaces, "diffusion operator Sff");
  LaplaceBeltramiFlat(surfaces, "diffusion operator Scc");
  LaplaceBeltramiFlat(surfaces, "diffusion operator");

  // internal surface(s)
  surfaces[0] = "Middle z-surface";
  LaplaceBeltramiFlat(surfaces, "diffusion operator");
  surfaces.push_back("Middle y-surface");
  LaplaceBeltramiFlat(surfaces, "diffusion operator");
  surfaces.push_back("Middle x-surface");
  LaplaceBeltramiFlat(surfaces, "diffusion operator");

  // half surfaces
  surfaces.resize(1);
  surfaces[0] = "Half z-surface";
  LaplaceBeltramiFlat(surfaces, "diffusion operator");
}


