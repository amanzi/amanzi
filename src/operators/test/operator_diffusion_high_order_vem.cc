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
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "SurfaceCoordinateSystem.hh"

// Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "VEM_Diffusion_HighOrder.hh"


/* *****************************************************************
* Exactness test for high-order Raviart-Thomas element.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_HIGH_ORDER_RAVIART_THOMAS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::WhetStone;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, high-order Raviart-Thomas" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  Teuchos::RCP<GeometricModel> gm;
  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, true, true);

  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create boundary data (no mixed bc)
  ParameterList op_list = plist.sublist("PK operator")
                               .sublist("diffusion operator Raviart-Thomas");
  int order = op_list.sublist("schema").get<int>("method order");
  int nd = PolynomialSpaceDimension(2, order);

  Point xv(3), x0(3), x1(3);
  AmanziMesh::Entity_ID_List nodes;

  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, DOF_Type::VECTOR));
  std::vector<int>& bc_model_f = bc_f->bc_model();
  std::vector<std::vector<double> >& bc_value_f = bc_f->bc_value_vector(nd);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {

      bc_model_f[f] = OPERATOR_BC_DIRICHLET;
      bc_value_f[f][0] = 0.0;
      bc_value_f[f][1] = 0.0;
      bc_value_f[f][2] = 0.0;
    }
  }

  // create diffusion operator 
  Teuchos::RCP<PDE_Abstract> op = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
  op->AddBCs(bc_f, bc_f);
  
  // populate the diffusion operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->Init();
std::cout << "HERE" << std::endl;
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
std::cout << "HERE" << std::endl;
global_op->rhs()->Print(std::cout);
  op->ApplyBCs(true, true, true);

  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditioner using the base operator class
  ParameterList slist;
  // slist.set<std::string>("preconditioner type", "diagonal");
  slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // solve the problem
  ParameterList lop_list = plist.sublist("solvers")
                                .sublist("AztecOO CG").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  solver.ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (pcg): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;
  }
exit(0);
}

