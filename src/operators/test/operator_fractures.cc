/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
#include "Tensor.hh"

// Amanzi::Operators
#include "Analytic02.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"


/* *****************************************************************
 * TBW.
 * **************************************************************** */
void
RunTest(int icase, double gravity)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "\nTest: Darcy flow in fractures, gravity=" << gravity
              << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_fractures.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> surfmesh;

  if (icase == 0) {
    RCP<const Mesh> mesh =
      meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

    // extract fractures mesh
    std::vector<std::string> setnames;
    setnames.push_back("fracture 1");
    setnames.push_back("fracture 2");
    surfmesh = meshfactory.create(mesh, setnames, AmanziMesh::FACE);
  } else {
    surfmesh = meshfactory.create("test/fractures.exo");
  }

  // modify diffusion coefficient
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned =
    surfmesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost =
    surfmesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_owned; c++) {
    WhetStone::Tensor Kc(2, 1);
    Kc(0, 0) = 1.0;
    K->push_back(Kc);
  }

  AmanziGeometry::Point v(1.0, 2.0, 3.0);
  Analytic02 ana(surfmesh, v, gravity);

  // create boundary data (no mixed bc)
  Teuchos::RCP<BCs> bc =
    Teuchos::rcp(new BCs(surfmesh, AmanziMesh::FACE, DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = surfmesh->face_centroid(f);
    if (fabs(xf[2] - 0.5) < 1e-6 &&
        (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
         fabs(xf[1] - 1.0) < 1e-6)) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    } else if (fabs(xf[1] - 0.5) < 1e-6 &&
               (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
                fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6)) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }

  // create solution
  Teuchos::RCP<CompositeVectorSpace> cvs =
    Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(surfmesh)->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1)->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector solution(*cvs);
  solution.putScalar(0.0);

  // create diffusion operator
  double rho(1.0);
  AmanziGeometry::Point g(0.0, 0.0, -gravity);
  Teuchos::ParameterList olist =
    plist.sublist("PK operator").sublist("diffusion operator");
  olist.set<bool>("gravity", (gravity > 0.0));

  Operators::PDE_DiffusionFactory opfactory;
  Teuchos::RCP<Operators::PDE_Diffusion> op =
    opfactory.Create(olist, surfmesh, bc, rho, g);
  op->SetBCs(bc, bc);

  Teuchos::RCP<Operator> global_op = op->global_operator();
  global_op->Init();

  // populate diffusion operator
  op->Setup(K, Teuchos::null, Teuchos::null);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);

  // apply BCs and assemble
  op->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  // create preconditoner
  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // solve the problem
  ParameterList lop_list =
    plist.sublist("solvers").sublist("PCG").sublist("pcg parameters");
  AmanziSolvers::
    LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector rhs = *global_op->rhs();
  int ierr = solver.applyInverse(rhs, solution);

  double a;
  a = rhs.norm2();
  if (getRank == 0) {
    std::cout << "pressure solver (" << solver.name()
              << "): ||r||=" << solver.residual()
              << " itr=" << solver.num_itrs() << "  ||f||=" << a
              << " code=" << solver.returned_code() << std::endl;
  }

  // calculate error in cell-centered data
  double pnorm, l2_err, inf_err;
  Epetra_MultiVector& p = *solution.ViewComponent("cell");

  ana.ComputeCellError(p, 0.0, pnorm, l2_err, inf_err);
  CHECK(l2_err < 1e-12);

  if (getRank == 0) {
    l2_err /= pnorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f\n", l2_err, inf_err);
  }

  // remove gravity to check symmetry
  if (gravity > 0.0) {
    for (int c = 0; c < ncells_owned; c++) {
      const Point& xc = surfmesh->cell_centroid(c);
      p[0][c] -= rho * g[2] * xc[2];
    }
  }

  if (getRank == 0) {
    GMV::open_data_file(*surfmesh, (std::string) "operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    GMV::close_data_file();
  }
}


TEST(FRACTURES_EXTRACTION)
{
  RunTest(0, 0.0);
}

TEST(FRACTURES_INPUT_EXODUS_FILE_GRAVITY)
{
  RunTest(1, 2.0);
}
