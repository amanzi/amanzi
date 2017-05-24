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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mfd3d_diffusion.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "DiffusionMFD.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"


/* *****************************************************************
* TBW.
* **************************************************************** */
void RunTest() {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Farcy flow in fractures" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_fractures.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, gm);
  // RCP<const Mesh> mesh = meshfactory("test/fractures.exo", gm);
  RCP<const Mesh_MSTK> mesh_mstk = rcp_static_cast<const Mesh_MSTK>(mesh);

  // extract surface mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  RCP<Mesh> surfmesh = Teuchos::rcp(new Mesh_MSTK(*mesh_mstk, setnames, AmanziMesh::FACE));

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  // create boundary data (no mixed bc)
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = surfmesh->face_centroid(f);
    if (fabs(xf[2] - 0.5) < 1e-6 && fabs(xf[0]) < 1e-6) { 
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = 1.0;
    }
    else if (fabs(xf[2] - 0.5) < 1e-6 && fabs(xf[0] - 1.0) < 1e-6) { 
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = 0.0;
    }
  }

  // create solution 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh)->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1)->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);

  // create diffusion operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
  DiffusionMFD op(olist, surfmesh);
  op.SetBCs(bc, bc);

  Teuchos::RCP<Operator> global_op = op.global_operator();
  global_op->Init();

  // populate diffusion operator
  op.Setup(K, Teuchos::null, Teuchos::null);
  op.UpdateMatrices(Teuchos::null, Teuchos::null);

  // apply BCs and assemble
  op.ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();
    
  // create preconditoner
  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers").sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  int ierr = solver.ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    double a;
    rhs.Norm2(&a);
    std::cout << "pressure solver (" << solver.name() 
              << "): ||r||=" << solver.residual() << " itr=" << solver.num_itrs()
              << "  ||f||=" << a 
              << " code=" << solver.returned_code() << std::endl;
  }
 
  if (MyPID == 0) {
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*surfmesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}


TEST(FRACTURES_EXTRACTION) {
  RunTest();
}


