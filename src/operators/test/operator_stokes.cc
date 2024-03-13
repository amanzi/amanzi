/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
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
#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "PDE_Accumulation.hh"
#include "PDE_Abstract.hh"
#include "PDE_Elasticity.hh"
#include "TreeOperator.hh"
#include "InverseFactory.hh"

#include "AnalyticElasticity02.hh"
#include "Verification.hh"


/* *****************************************************************
* Stokes model: exactness test.
***************************************************************** */
SUITE(OPERATOR_STOKES)
{
  TEST(OPERATOR_STOKES_EXACTNESS_FORWARD)
  {
    using namespace Amanzi;
    using namespace Amanzi::AmanziMesh;
    using namespace Amanzi::AmanziGeometry;
    using namespace Amanzi::Operators;

    auto comm = Amanzi::getDefaultComm();
    int MyPID = comm->MyPID();
    if (MyPID == 0) std::cout << "\nTest: 2D Stokes: exactness test" << std::endl;

    // read parameter list
    std::string xmlFileName = "test/operator_stokes.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();

    // create a simple rectangular mesh
    MeshFactory meshfactory(comm);
    meshfactory.set_preference(Preference({ Framework::MSTK }));
    Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 12, 14);

    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    int nfaces_wghost =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
    int nnodes_wghost =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

    // populate diffusion coefficient using an analytic solution
    AnalyticElasticity02 ana(mesh);

    Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
      Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->getCellCentroid(c);
      const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
      K->push_back(Kc);
    }

    // populate boundary data for Bernardi-Raugel-type element
    // -- Dirichlet condition on faces for the normal velocity component
    Teuchos::RCP<BCs> bcf =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bcf_model = bcf->bc_model();
    std::vector<double>& bcf_value = bcf->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
          fabs(xf[1] - 1.0) < 1e-6) {
        const Point& normal = mesh->getFaceNormal(f);
        double area = mesh->getFaceArea(f);

        bcf_model[f] = OPERATOR_BC_DIRICHLET;
        bcf_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
      }
    }

    // -- Dirichlet condition at nodes for the normal velocity component
    Point xv(2);
    Teuchos::RCP<BCs> bcv =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
    std::vector<int>& bcv_model = bcv->bc_model();
    std::vector<Point>& bcv_value = bcv->bc_value_point();

    for (int v = 0; v < nnodes_wghost; ++v) {
      xv = mesh->getNodeCoordinate(v);

      if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1]) < 1e-6 ||
          fabs(xv[1] - 1.0) < 1e-6) {
        bcv_model[v] = OPERATOR_BC_DIRICHLET;
        bcv_value[v] = ana.velocity_exact(xv, 0.0);
      }
    }

    // create a discrete PDE
    // -- create an elasticity operator using in particular schema data. In the future,
    // -- we could take definition of DOFs data from WhetStone using the provided
    // -- discretization method.
    Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator");
    Teuchos::RCP<PDE_Elasticity> op00 = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
    op00->SetTensorCoefficient(K);
    // -- add boundary conditions to the discrete PDE: Lame equation
    op00->SetBCs(bcf, bcf);
    op00->AddBCs(bcv, bcv);

    // -- create a divergence operator (incompressibility equation).
    //    Note that it uses two different schemas for DOFs for velocity and pressure.
    op_list = plist.sublist("PK operator").sublist("divergence operator");
    Teuchos::RCP<PDE_Abstract> op10 = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
    op10->SetBCs(bcf, bcf);
    op10->AddBCs(bcv, bcv);

    // -- create a gradient operator as transpose of divergence
    op_list = plist.sublist("PK operator").sublist("gradient operator");
    Teuchos::RCP<PDE_Abstract> op01 = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
    op01->SetBCs(bcf, bcf);
    op01->AddBCs(bcv, bcv);

    // create a tree operator for the discrete Stokes PDE. It is a 2x2 block operator
    // composed of Lame operator (blok 00), divergence operator (blok 10) and traspose
    // of the divergence operator (block 01).
    Teuchos::RCP<Operator> global00 = op00->global_operator();
    Teuchos::RCP<Operator> global10 = op10->global_operator();
    Teuchos::RCP<Operator> global01 = op01->global_operator();

    const CompositeVectorSpace& cvs = global00->DomainMap();
    Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
    tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs))));
    tvs->PushBack(
      Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op10->global_operator()->RangeMap()))));

    Teuchos::RCP<TreeOperator> op = Teuchos::rcp(new Operators::TreeOperator(tvs));
    op->set_operator_block(0, 0, global00);
    op->set_operator_block(1, 0, global10);
    op->set_operator_block(0, 1, global01);

    // create source term representing external forces.
    CompositeVector source(cvs);
    Epetra_MultiVector& src = *source.ViewComponent("node");

    for (int v = 0; v < nnodes; v++) {
      xv = mesh->getNodeCoordinate(v);
      Point tmp(ana.source_exact(xv, 0.0));
      for (int k = 0; k < 2; ++k) src[k][v] = tmp[k];
    }

    // populate local matrices inside the elasticity block and assemble them
    // into a global matrix.
    op00->UpdateMatrices();
    global00->UpdateRHS(source, true);
    op00->ApplyBCs(true, true, true);

    // populate local matrices inside the divergence block. Since we will use
    // a matrix-free matvec inside an iterative solver, there is no need to
    // assemble a global matrix.
    op10->UpdateMatrices();
    op10->ApplyBCs(false, true, false);

    op01->UpdateMatrices();
    op01->ApplyBCs(true, false, false);

    // -- copy right-hand sides inside two operators to the global rhs.
    TreeVector rhs(*tvs);
    *rhs.SubVector(0)->Data() = *global00->rhs();
    *rhs.SubVector(1)->Data() = *global10->rhs();

    // create and set solution variables: velocity and pressure.
    TreeVector solution(*tvs);
    ana.VectorNodeSolution(*solution.SubVector(0)->Data(), 0.0);
    ana.ScalarFaceSolution(*solution.SubVector(0)->Data(), 0.0);
    ana.ScalarCellSolution(*solution.SubVector(1)->Data(), 0.0);

    TreeVector res(*tvs);
    op->Apply(solution, res);
    res.Update(-1.0, rhs, 1.0);

    // Post-processing
    // -- compute velocity resor
    double ul2_res, uinf_res, unorm;
    res.SubVector(0)->Data()->Norm2(&ul2_res);
    res.SubVector(0)->Data()->NormInf(&uinf_res);
    solution.SubVector(0)->Data()->Norm2(&unorm);

    // -- compute pressure error
    double pl2_res, pinf_res;
    res.SubVector(1)->Data()->Norm2(&pl2_res);
    res.SubVector(1)->Data()->NormInf(&pinf_res);

    if (MyPID == 0) {
      ul2_res /= unorm;
      printf("L2(u)=%12.8g  Inf(u)=%12.8g  L2(p)=%12.8g  Inf(p)=%12.8g\n",
             ul2_res,
             uinf_res,
             pl2_res,
             pinf_res);

      CHECK(ul2_res < 0.01);
      CHECK(pl2_res < 0.0001);
    }
  }


  /* *****************************************************************
* Stokes model: exactness test.
***************************************************************** */
  TEST(OPERATOR_STOKES_EXACTNESS_INVERSE)
  {
    using namespace Amanzi;
    using namespace Amanzi::AmanziMesh;
    using namespace Amanzi::AmanziGeometry;
    using namespace Amanzi::Operators;

    auto comm = Amanzi::getDefaultComm();
    int MyPID = comm->MyPID();
    if (MyPID == 0) std::cout << "\nTest: 2D Stokes: exactness test" << std::endl;

    // read parameter list
    std::string xmlFileName = "test/operator_stokes.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();

    // create a simple rectangular mesh
    MeshFactory meshfactory(comm);
    meshfactory.set_preference(Preference({ Framework::MSTK }));
    Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 16, 20);

    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    int nfaces_wghost =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
    int nnodes_wghost =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

    // populate diffusion coefficient using an analytic solution
    AnalyticElasticity02 ana(mesh);

    Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
      Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->getCellCentroid(c);
      const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
      K->push_back(Kc);
    }

    // populate boundary data for Bernardi-Raugel-type element
    // -- Dirichlet condition on faces for the normal velocity component
    Teuchos::RCP<BCs> bcf =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    std::vector<int>& bcf_model = bcf->bc_model();
    std::vector<double>& bcf_value = bcf->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      const Point& normal = mesh->getFaceNormal(f);
      double area = mesh->getFaceArea(f);

      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
          fabs(xf[1] - 1.0) < 1e-6) {
        bcf_model[f] = OPERATOR_BC_DIRICHLET;
        bcf_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
      }
    }

    // -- Dirichlet condition at nodes for the normal velocity component
    Point xv(2);
    Teuchos::RCP<BCs> bcv =
      Teuchos::rcp(new BCs(mesh, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
    std::vector<int>& bcv_model = bcv->bc_model();
    std::vector<Point>& bcv_value = bcv->bc_value_point();

    for (int v = 0; v < nnodes_wghost; ++v) {
      xv = mesh->getNodeCoordinate(v);

      if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1]) < 1e-6 ||
          fabs(xv[1] - 1.0) < 1e-6) {
        bcv_model[v] = OPERATOR_BC_DIRICHLET;
        bcv_value[v] = ana.velocity_exact(xv, 0.0);
      }
    }

    // create a discrete PDE
    // -- create an elasticity operator using in particular schema data. In the future,
    // -- we could take definition of DOFs data from WhetStone using the provided
    // -- discretization method.
    Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist("elasticity operator");
    Teuchos::RCP<PDE_Elasticity> op00 = Teuchos::rcp(new PDE_Elasticity(op_list, mesh));
    op00->SetTensorCoefficient(K);
    // -- add boundary conditions to the discrete PDE: Lame equation
    op00->SetBCs(bcf, bcf);
    op00->AddBCs(bcv, bcv);

    // -- create a divergence operator (incompressibility equation).
    //    Note that it uses two different schemas for DOFs for velocity and pressure.
    op_list = plist.sublist("PK operator").sublist("divergence operator");
    Teuchos::RCP<PDE_Abstract> op10 = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
    op10->SetBCs(bcf, bcf);
    op10->AddBCs(bcv, bcv);

    // -- create a gradient operator as transpose of divergence
    op_list = plist.sublist("PK operator").sublist("gradient operator");
    Teuchos::RCP<PDE_Abstract> op01 = Teuchos::rcp(new PDE_Abstract(op_list, mesh));
    op01->SetBCs(bcf, bcf);
    op01->AddBCs(bcv, bcv);

    // create identity type operator: pressure block for preconditioner
    Teuchos::RCP<PDE_Accumulation> pc11 =
      Teuchos::rcp(new PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh));
    Teuchos::RCP<Operator> global11 = pc11->global_operator();

    // create a tree operator for the discrete Stokes PDE. It is a 2x2 block operator
    // composed of Lame operator (blok 00), divergence operator (blok 10) and traspose
    // of the divergence operator (block 01).
    Teuchos::RCP<Operator> global00 = op00->global_operator();
    Teuchos::RCP<Operator> global10 = op10->global_operator();
    Teuchos::RCP<Operator> global01 = op01->global_operator();

    const CompositeVectorSpace& cvs = global00->DomainMap();
    Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
    tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(cvs))));
    tvs->PushBack(
      Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op10->global_operator()->RangeMap()))));

    Teuchos::RCP<TreeOperator> op = Teuchos::rcp(new Operators::TreeOperator(tvs));
    op->set_operator_block(0, 0, global00);
    op->set_operator_block(1, 0, global10);
    op->set_operator_block(0, 1, global01);

    // create and initialize state variables: velocity and pressure.
    TreeVector solution(*tvs);
    solution.PutScalar(0.0);

    // create source term representing external forces.
    CompositeVector source(cvs);
    Epetra_MultiVector& src = *source.ViewComponent("node");

    for (int v = 0; v < nnodes; v++) {
      xv = mesh->getNodeCoordinate(v);
      Point tmp(ana.source_exact(xv, 0.0));
      for (int k = 0; k < 2; ++k) src[k][v] = tmp[k];
    }

    // populate local matrices inside the elasticity block and assemble them
    // into a global matrix.
    op00->UpdateMatrices();
    global00->UpdateRHS(source, true);
    op00->ApplyBCs(true, true, true);

    // populate local matrices inside the divergence block. Since we will use
    // a matrix-free matvec inside an iterative solver, there is no need to
    // assemble a global matrix.
    op10->UpdateMatrices();
    op10->ApplyBCs(false, true, false);

    op01->UpdateMatrices();
    op01->ApplyBCs(true, false, false);

    // populate local matrices in the pressure block (for preconditioner)
    CompositeVector vol(global11->DomainMap());
    vol.PutScalar(1.0);
    pc11->AddAccumulationTerm(vol, 1.0, "cell");

    // create a block-diagonal preconditoner, identity is the default one.
    // The first block will reuse the assembled matrix to build a multigrid
    // solver. The second block is simply a diagonal matrix. The off-diagonal
    // blocks are empty and require no setup.
    Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
    global00->set_inverse_parameters(slist);

    slist = plist.sublist("preconditioners").sublist("Diagonal");
    global11->set_inverse_parameters(slist);

    Teuchos::RCP<TreeOperator> pc = Teuchos::rcp(new Operators::TreeOperator(tvs));
    pc->set_operator_block(0, 0, op00->global_operator());
    pc->set_operator_block(1, 1, pc11->global_operator());

    Teuchos::ParameterList solver_list;
    solver_list.set("preconditioning method", "block diagonal");
    pc->set_inverse_parameters(solver_list);
    pc->InitializeInverse();
    pc->ComputeInverse();

    // Assemble global matrix for tesing purposes
    // // Test SPD properties of the matrix and preconditioner.
    // VerificationTV ver1(op), ver2(pc);
    // ver1.CheckMatrixSPD(true, false, false);
    // ver1.CheckMatrixSPD(true, false, true);
    // ver2.CheckPreconditionerSPD();

    // now set up a full solver
    //
    // NOTE: this cannot be done in any TreeOperator, as we are using _different
    // structure_ for the forward operator and preconditioner.  Instead we
    // directly call the factory.
    solver_list.setParameters(plist.sublist("solvers").sublist("GMRES"));
    auto solver = AmanziSolvers::createIterativeMethod(solver_list, op, pc);
    solver->InitializeInverse();
    solver->ComputeInverse();

    // -- copy right-hand sides inside two operators to the global rhs.
    TreeVector rhs(*tvs);
    *rhs.SubVector(0)->Data() = *global00->rhs();
    *rhs.SubVector(1)->Data() = *global10->rhs();

    // -- execute GMRES solver
    solver->ApplyInverse(rhs, solution);

    // solver->AssembleMatrix();
    // ver1.CheckResidual(solution, rhs, 1.0e-12);

    if (MyPID == 0) {
      std::cout << "elasticity solver (gmres): ||r||=" << solver->residual()
                << " itr=" << solver->num_itrs() << " code=" << solver->returned_code()
                << std::endl;
    }

    // Post-processing
    // -- compute velocity error
    double unorm, ul2_err, uinf_err;
    ana.VectorNodeError(*solution.SubVector(0)->Data(), 0.0, unorm, ul2_err, uinf_err);

    // -- compute pressure error
    double pnorm, pl2_err, pinf_err;
    ana.ScalarCellError(*solution.SubVector(1)->Data(), 0.0, pnorm, pl2_err, pinf_err);

    if (MyPID == 0) {
      ul2_err /= unorm;
      printf("L2(u)=%12.8g  Inf(u)=%12.8g  L2(p)=%12.8g  Inf(p)=%12.8g  itr=%3d\n",
             ul2_err,
             uinf_err,
             pl2_err,
             pinf_err,
             solver->num_itrs());

      CHECK(ul2_err < 0.01);
      CHECK(pl2_err < 0.05);
      CHECK(solver->num_itrs() < 60);
    }

    if (MyPID == 0) {
      const Epetra_MultiVector& u = *solution.SubVector(0)->Data()->ViewComponent("node");
      GMV::open_data_file(*mesh, (std::string) "operators.gmv");
      GMV::start_data();
      GMV::write_node_data(u, 0, "velocity_x");
      GMV::write_node_data(u, 1, "velocity_y");
      GMV::close_data_file();
    }
  }

} // SUITE
