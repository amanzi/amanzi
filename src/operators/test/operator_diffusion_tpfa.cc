/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

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
#include "LinearOperatorFactory.hh"
#include "mfd3d_diffusion.hh"
#include "tensor.hh"

// Operators
#include "Analytic04.hh"
#include "HeatConduction.hh"

#include "DiffusionTPFA.hh"
#include "OperatorDefs.hh"
#include "OperatorSource.hh"
#include "UpwindSecondOrder.hh"
#include "UpwindStandard.hh"


/* *****************************************************************
* This tests diffusion solvers with zero coefficients
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_TPFA_ZEROCOEF)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0)
    std::cout << "\nTest: 2D elliptic solver, TPFA with zero permeability" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_strip.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory(-4.0, 0.0, 4.0, 1.0, 30, 1, gm);

  // model
  Analytic04 ana(mesh);

  // vector spaces
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);

  Teuchos::RCP<CompositeVectorSpace> cell_space = Teuchos::rcp(new CompositeVectorSpace());
  cell_space->SetMesh(mesh)->SetComponent("cell", CELL, 1)->SetGhosted();

  Teuchos::RCP<CompositeVectorSpace> face_space = Teuchos::rcp(new CompositeVectorSpace());
  face_space->SetMesh(mesh)->SetComponent("face", FACE, 1)->SetGhosted();

  // create source
  CompositeVector source(*cell_space);
  Epetra_MultiVector& src = *source.ViewComponent("cell", false);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->getCellCentroid(c);
    double volume = mesh->getCellVolume(c);
    src[0][c] = ana.source_exact(xc, 0.0) * volume;
  }

  Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(cell_space, 0));
  op1->UpdateMatrices(source);

  // modify diffusion coefficient
  // -- since rho=mu=1.0, we do not need to scale the tensor coefficient.
  std::vector<WhetStone::Tensor> K;

  for (int c = 0; c != ncells; ++c) {
    const Point& xc = mesh->getCellCentroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K.push_back(Kc);
  }

  // -- scalar part
  Teuchos::RCP<CompositeVector> coef = Teuchos::rcp(new CompositeVector(*face_space));
  {
    Epetra_MultiVector& coef_faces = *coef->ViewComponent("face", false);

    for (int f = 0; f != nfaces; ++f) {
      const Point& xf = mesh->getFaceCentroid(f);
      coef_faces[0][f] = ana.ScalarCoefficient(xf, 0.0);
    }
  }
  coef->ScatterMasterToGhosted("face");

  // create boundary data
  Point xv(2);
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);
  std::vector<double> bc_mixed;

  AmanziMesh::Entity_ID_List left;
  mesh->getSetEntities("Left side", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &left);
  for (int f = 0; f != left.size(); ++f) {
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    mesh->getFaceCentroid(f, &xv);
    bc_value[f] = ana.pressure_exact(xv, 0.0);
  }

  AmanziMesh::Entity_ID_List right;
  mesh->getSetEntities("Right side", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &right);
  for (int f = 0; f != right.size(); ++f) {
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    mesh->getFaceCentroid(f, &xv);
    bc_value[f] = ana.pressure_exact(xv, 0.0);
  }

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator
  ParameterList op_list = plist.sublist("PK operator").sublist("diffusion operator");
  Point g(0.0, 0.0);

  Teuchos::RCP<DiffusionTPFA> op2 = Teuchos::rcp(new DiffusionTPFA(*op1, op_list, bc));
  op2->SetUpwind(0);
  op2->SetGravity(g);

  const CompositeVectorSpace& cvs = op2->DomainMap();

  // populate the diffusion operator
  op2->Setup(K, coef, Teuchos::null);
  op2->UpdateMatrices(Teuchos::null, Teuchos::null);
  op2->ApplyBCs();
  op2->SymbolicAssembleMatrix(op2->schema_prec_dofs());
  op2->AssembleMatrix(op2->schema_prec_dofs());

  // create solution vector
  CompositeVector solution(*cell_space);
  {
    Epetra_MultiVector& p_cell = *solution.ViewComponent("cell");
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->getCellCentroid(c);
      p_cell[0][c] = ana.pressure_exact(xc, 0.0);
    }
  }

  // check residual
  CompositeVector residual(*cell_space);
  // int ierr = op1->ComputeNegativeResidual(solution, residual);
  solution.Print(std::cout);

  op2->rhs()->Print(std::cout);
  int ierr = op2->Apply(solution, residual);

  residual.Print(std::cout);
  residual.Update(1.0, *op2->rhs(), -1.0);
  residual.Print(std::cout);
  CHECK(!ierr);

  double res_norm(0.0);
  ierr |= residual.Norm2(&res_norm);
  CHECK(!ierr);
  CHECK_CLOSE(0.0, res_norm, 1.0e-8);
}
