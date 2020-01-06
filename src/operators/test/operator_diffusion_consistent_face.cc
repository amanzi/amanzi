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
#include "EpetraExt_RowMatrixOut.h"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "Tensor.hh"

// Operators
#include "Analytic02.hh"

#include "OperatorDefs.hh"
#include "PDE_DiffusionMFD.hh"


/* *****************************************************************
 * This does a check that the UpdateConsistentFace method results in face
 * values that satisfy the linear equation.
 * **************************************************************** */
TEST(OPERATOR_DIFFUSION_MIXED) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
                            << " test for mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_consistent_face.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  //RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 1);
  RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the diffusion coefficient.
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic02 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();
  std::vector<double>& bc_mixed = bc->bc_mixed();

  for (int f = 0; f < nfaces_wghost; f++) {
    bool flag;
    const Point& xf = mesh->face_centroid(f);
    double area = mesh->face_area(f);
    Point normal = ana.face_normal_exterior(f, &flag);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;  // We assume exterior normal.
    } else if (fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;  // We assume exterior normal.

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    } else if (fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create diffusion operator 
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator mixed");
  auto op = Teuchos::rcp(new PDE_DiffusionMFD(op_list, mesh));
  op->Init();
  op->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op->global_operator()->DomainMap();

  // set up the diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // get and assmeble the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();
  op->ApplyBCs(true, true, true);

  // put in a random set of cell values
  CompositeVector x(cvs);
  x.Random();
  x.ViewComponent("face")->PutScalar(0.);

  op->UpdateConsistentFaces(x);

  // // dump the schur complement
  // std::stringstream filename_s2;
  // filename_s2 << "consist_face_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *op->consistent_face_operator()->A());

  // ensure that (y - A * x) on faces is zero
  CompositeVector res(cvs);
  global_op->ComputeNegativeResidual(x,res);

  double norm;
  res.ViewComponent("face",false)->NormInf(&norm);
  CHECK_CLOSE(0.0, norm, 1.e-8);
}
