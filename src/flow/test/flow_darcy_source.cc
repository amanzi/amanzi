/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
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

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"
#include "LinearOperatorFactory.hh"

#include "tensor.hh"
#include "mfd3d_diffusion.hh"

#include "State.hh"
#include "Flow_State.hh"
#include "Darcy_PK.hh"
#include "Stiffness_MFD.hh"


double Kxx(const Amanzi::AmanziGeometry::Point& p, double t) {
  double x = p[0];
  double y = p[1]; 
  return (x + 1) * (x + 1) + y * y;
}
double Kyy(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];
  return (x + 1) * (x + 1);
}
double Kxy(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];
  return -x * y;
}
double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];
  double xy = x * y;
  return x * xy * xy + x * sin(2 * M_PI * xy) * sin(2 * M_PI * y);
}
Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];

  double t01, t02, t03, t12, t13, t04, t05, t06; 
  double px, py;

  t01 = x*x*y;
  t02 = sin(2*M_PI*x*y);
  t03 = sin(2*M_PI*y);

  t12 = cos(2*M_PI*x*y);
  t13 = cos(2*M_PI*y);

  px = 3*y*t01 + t03*(t02 + 2*M_PI*y*x*t12);
  py = 2*x*t01 + x*2*M_PI*(x*t12*t03 + t02*t13);

  t04 = Kxx(p, t);
  t05 = Kxy(p, t);
  t06 = Kyy(p, t);

  Amanzi::AmanziGeometry::Point v(2);
  v[0] = -(t04 * px + t05 * py);
  v[1] = -(t05 * px + t06 * py);
  return v;
}
Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];

  double t01, t02, t03, t12, t13, t04, t05, t06; 
  double px, py;

  t01 = x*x*y;
  t02 = sin(2*M_PI*x*y);
  t03 = sin(2*M_PI*y);

  t12 = cos(2*M_PI*x*y);
  t13 = cos(2*M_PI*y);

  px = 3*y*t01 + t03*(t02 + 2*M_PI*y*x*t12);
  py = 2*x*t01 + x*2*M_PI*(x*t12*t03 + t02*t13);

  Amanzi::AmanziGeometry::Point v(2);
  v[0] = px;
  v[1] = py;
  return v;
}
double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];

  double t01, t02, t03, t12, t13;
  double px, py, pxx, pxy, pyy;
  double t04, t05, t06, tx4, ty4, tx5, ty5, tx6;

  t01 = x*x*y;
  t02 = sin(2*M_PI*x*y);
  t03 = sin(2*M_PI*y);

  t12 = cos(2*M_PI*x*y);
  t13 = cos(2*M_PI*y);

  px = 3*y*t01 + t03*(t02 + 2*M_PI*y*x*t12);
  py = 2*x*t01 + x*2*M_PI*(x*t12*t03 + t02*t13);

  pxx = 6*x*y*y + 4*M_PI*t03*(y*t12 - M_PI*y*y*x*t02); 
  pxy = 6*x*x*y + 2*M_PI*(t13*t02 + 2*x*t03*t12 + x*y*2*M_PI*(t13*t12-x*t03*t02));
  pyy = 2*x*x*x + x*4*M_PI*M_PI*(-x*x*t02*t03 + 2*x*t12*t13 - t02*t03);

  t04 = Kxx(p, t);
  t05 = Kxy(p, t);
  t06 = Kyy(p, t);

  tx4 = 2*(x+1);  // d/dx (Kxx)
  ty4 = 2*y;      // d/dy (Kxx)

  tx5 = -y;  // d/dx (Kxy)  
  ty5 = -x;  // d/dy (Kxy)
  
  tx6 = 2*(x+1);  // d/dy (Kxy)
  return -(tx4 + ty5)*px - tx5*py - t04*pxx - 2*t05*pxy - t06*pyy;
}

/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. It is most sensitive to error in location of problem
* coefficients.
* **************************************************************** */
TEST(FLOW_DARCY_SOURCE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "Test: 2D steady-state elliptic solver, mixed discretization" << endl;

  /* read parameter list */
  string xmlFileName = "test/flow_darcy_source.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create an SIMPLE mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 200, 200, gm);
  RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Darcy_PK* DPK = new Darcy_PK(plist, S);
  S->Setup();
  S->InitializeFields();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Epetra_MultiVector& perm = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  
  for (int c = 0; c < perm.MyLength(); c++) {
    perm[0][c] = 2.0;
    perm[1][c] = 1.0;
  }

  *S->GetScalarData("fluid_density", passwd) = 1.0;
  *S->GetScalarData("fluid_viscosity", passwd) = 1.0;
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", passwd);
  gravity[1] = 0.0;
  S->GetFieldData("porosity", passwd)->PutScalar(1.0);
 
  /* initialize the Darcy process kernel */
  DPK->InitPK();
  DPK->InitSteadyState(0.0, 1e-8);

  /* reset problem coefficients */
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  std::vector<WhetStone::Tensor>& K = DPK->get_K();

  for (int c = 0; c < K.size(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    K[c](0, 0) = Kxx(xc, 0.0);
    K[c](1, 1) = Kyy(xc, 0.0);
    K[c](0, 1) = Kxy(xc, 0.0);
    K[c](1, 0) = Kxy(xc, 0.0);
  }

  /* update Dirichlet boundary data */
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  std::vector<bc_tuple>& bc_values = DPK->get_bc_values();

  for (int f = 0; f < nfaces; f++) {
    if (bc_values[f][0] > 0.0) {
      const Point& xf = mesh->face_centroid(f);
      bc_values[f][0] = pressure_exact(xf, 0.0);
    }
  }

  for (int n = 0; n < 240; n+=50) {
    // recalculate mass matrices
    double factor = pow(10.0, (double)(n - 50) / 100.0);
    Matrix_MFD* matrix = dynamic_cast<Matrix_MFD*>(&*DPK->matrix());
    matrix->CreateMassMatrices_ScaledStability(AmanziFlow::FLOW_MFD3D_POLYHEDRA, factor, K);

    // create a problem
    DPK->AssembleMatrixMFD();

    // update right-hand side
    Epetra_MultiVector& b = *DPK->matrix()->rhs()->ViewComponent("cell");
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      double volume = mesh->cell_volume(c);
      b[0][c] += source_exact(xc, 0.0) * volume;
    }

    // steady-state solution
    CompositeVector& rhs = *DPK->matrix()->rhs();
    CompositeVector& solution = *DPK->get_solution();
    DPK->SolveFullySaturatedProblem(0.0, rhs, solution);
    DPK->CommitState(S);

    // calculate errors
    Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
    Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", true);

    double p_norm(0.0), p_error(0.0);
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      double tmp = pressure_exact(xc, 0.0);
      double volume = mesh->cell_volume(c);

      p_error += std::pow(tmp - p[0][c], 2.0) * volume;
      p_norm += std::pow(tmp, 2.0) * volume;
    }

    double flux_norm(0.0), flux_error(0.0);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      const Point& xf = mesh->face_centroid(f);
      const AmanziGeometry::Point& velocity = velocity_exact(xf, 0.0);
      double tmp = velocity * normal;
 
      flux_error += std::pow(tmp - flux[0][f], 2.0);
      flux_norm += std::pow(tmp, 2.0);
    }  

    p_error = pow(p_error / p_norm, 0.5);
    flux_error = pow(flux_error / flux_norm, 0.5);
    printf("scale = %10.6g  Err(p) = %10.6f  Err(flux) = %10.6g\n", factor, p_error, flux_error); 
    
    CHECK(p_error< 0.15 && flux_error < 0.15);
  }

  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(p, 0, "pressure");
  GMV::close_data_file();

  delete DPK;
}


/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. This is a prototype for future solvers.
* **************************************************************** */
TEST(FLOW_DARCY_NODAL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "\nTest: 2D steady-state elliptic solver, nodal discretization" << endl;

  // read parameter list
  string xmlFileName = "test/flow_darcy_source.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 80, 80, gm);
  RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = Kxx(xc, 0.0);
    Kc(1, 1) = Kyy(xc, 0.0);
    Kc(0, 1) = Kxy(xc, 0.0);
    Kc(1, 0) = Kxy(xc, 0.0);

    K.push_back(Kc);
  }

  // create boundary data
  Point xv(2);
  std::vector<int> bc_model(nnodes);
  std::vector<double> bc_values(nnodes);

  for (int v = 0; v < nnodes; v++) {
    mesh->node_get_coordinates(v, &xv);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model[v] = FLOW_BC_FACE_PRESSURE;
      bc_values[v] = pressure_exact(xv, 0.0);
    }
  }

  // create matrix and preconditioner
  Teuchos::RCP<Stiffness_MFD> matrix = Teuchos::rcp(new Stiffness_MFD(mesh));
  Teuchos::RCP<Stiffness_MFD> preconditioner = Teuchos::rcp(new Stiffness_MFD(mesh));

  matrix->SymbolicAssembleGlobalMatrices();
  preconditioner->SymbolicAssembleGlobalMatrices();

  ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners")
                             .get<Teuchos::ParameterList>("Hypre AMG")
                             .get<Teuchos::ParameterList>("boomer amg parameters");
  preconditioner->InitPreconditioner(FLOW_PRECONDITIONER_HYPRE_AMG, slist);

  for (int n = 0; n < 400; n+=110) {
    // double factor = pow(10.0, (double)(n - 50) / 100.0);
    double factor = pow(10.0, (double)(n - 150) / 100.0);

    // populate matrix
    matrix->CreateMFDstiffnessMatrices(K, factor);
    matrix->CreateMFDrhsVectors();
    matrix->ApplyBoundaryConditions(bc_model, bc_values);
    matrix->AssembleGlobalMatrices();

    // populate preconditioner
    preconditioner->CreateMFDstiffnessMatrices(K, factor);
    preconditioner->CreateMFDrhsVectors();
    preconditioner->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner->AssembleGlobalMatrices();
    preconditioner->UpdatePreconditioner();
  
    // update right-hand side
    AmanziMesh::Entity_ID_List nodes;
    Teuchos::RCP<Epetra_Vector> rhs = matrix->rhs();

    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      double volume = mesh->cell_volume(c);

      mesh->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int k = 0; k < nnodes; k++) {
        int v = nodes[k];
        (*rhs)[v] += source_exact(xc, 0.0) * volume / nnodes;
      }
    }

    // solve the problem
    ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
    AmanziSolvers::LinearOperatorFactory<Stiffness_MFD, Epetra_Vector, Epetra_BlockMap> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Stiffness_MFD, Epetra_Vector, Epetra_BlockMap> >
       solver = factory.Create("AztecOO CG", lop_list, matrix, preconditioner);

    Epetra_Vector solution(*rhs);
    solver->ApplyInverse(*rhs, solution);

    // calculate errors
    WhetStone::MFD3D_Diffusion mfd(mesh);
    AmanziGeometry::Point grad(2);

    double p_norm(0.0), p_error(0.0);
    double grad_norm(0.0), grad_error(0.0);

    for (int c = 0; c < ncells; c++) {
      double volume = mesh->cell_volume(c);

      mesh->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      std::vector<double> cell_solution(nnodes);

      for (int k = 0; k < nnodes; k++) {
        int v = nodes[k];
        cell_solution[k] = solution[v];

        mesh->node_get_coordinates(v, &xv);
        double tmp = pressure_exact(xv, 0.0);

        p_error += std::pow(tmp - solution[v], 2.0) * volume / nnodes;
        p_norm += std::pow(tmp, 2.0) * volume / nnodes;
      }

      const Point& xc = mesh->cell_centroid(c);
      const AmanziGeometry::Point& grad_exact = gradient_exact(xc, 0.0);
      mfd.RecoverGradient_StiffnessMatrix(c, cell_solution, grad);

      grad_error += L22(grad - grad_exact) * volume;
      grad_norm += L22(grad_exact) * volume;
    }

    p_error = pow(p_error / p_norm, 0.5);
    grad_error = pow(grad_error / grad_norm, 0.5);
    printf("scale = %10.6g  Err(p) = %10.6f  Err(grad) = %10.6g\n", factor, p_error, grad_error); 

    CHECK(p_error < 0.1 && grad_error < 0.15);
  }
}
