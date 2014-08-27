/*
  This is the operator component of the Amanzi code. 

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

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"
#include "OperatorSource.hh"


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
TEST(OPERATOR_MIXED_DIFFUSION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D steady-state elliptic solver, mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::hide_line_prefix = true;

  // create a mesh 
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 4, 4, gm);
  RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);
  // RCP<const Mesh> mesh = meshfactory("test/median255x256.exo", gm);

  // create diffusion coefficient
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    WhetStone::Tensor Kc(2, 2);
    const Point& xc = mesh->cell_centroid(c);

    Kc(0, 0) = Kxx(xc, 0.0);
    Kc(1, 1) = Kyy(xc, 0.0);
    Kc(0, 1) = Kxy(xc, 0.0);
    Kc(1, 0) = Kxy(xc, 0.0);

    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);

  // create boundary data
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  Point xv(2);
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_value[f] = pressure_exact(xf, 0.0);
      bc_model[f] = Operators::OPERATOR_BC_FACE_DIRICHLET;
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(bc_model, bc_value));

  // create diffusion operator 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  cvs->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector solution(*cvs), flux(*cvs);
  solution.PutScalar(0.0);

  // create source 
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);
  
  Epetra_MultiVector& src = *source.ViewComponent("cell");
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);
    src[0][c] += source_exact(xc, 0.0);
  }

  // MAIN LOOP
  for (int n = 0; n < 240; n+=50) {
    double factor = pow(10.0, (double)(n - 50) / 100.0);
    
    // create source operator 
    Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(cvs, 0));
    op1->Init();
    op1->UpdateMatrices(source);

    // populate the diffusion operator
    Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operators")
                                        .get<Teuchos::ParameterList>("mixed diffusion");
    Teuchos::RCP<OperatorDiffusion> op2 = Teuchos::rcp(new OperatorDiffusion(*op1, olist, bc));

    int schema_dofs = op2->schema_dofs();
    int schema_prec_dofs = op2->schema_prec_dofs();
    CHECK(schema_dofs == Operators::OPERATOR_SCHEMA_DOFS_FACE + Operators::OPERATOR_SCHEMA_DOFS_CELL);
    // CHECK(schema_prec_dofs == Operators::OPERATOR_SCHEMA_DOFS_FACE);

    op2->set_factor(factor);  // for developers only
    op2->InitOperator(K, Teuchos::null, Teuchos::null, rho, mu);
    op2->UpdateMatrices(Teuchos::null, Teuchos::null);
    op2->ApplyBCs();
    op2->SymbolicAssembleMatrix(schema_prec_dofs);
    op2->AssembleMatrix(schema_prec_dofs);
    
    ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
    op2->InitPreconditioner("Hypre AMG", slist);

    // solve the problem
    ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
    solution.PutScalar(0.0);
    AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
       solver = factory.Create("AztecOO CG", lop_list, op2);

    CompositeVector& rhs = *op2->rhs();
    int ierr = solver->ApplyInverse(rhs, solution);

    // calculate pressure errors
    Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

    double p_norm(0.0), p_error(0.0);
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      double tmp = pressure_exact(xc, 0.0);
      double volume = mesh->cell_volume(c);

      p_error += std::pow(tmp - p[0][c], 2.0) * volume;
      p_norm += std::pow(tmp, 2.0) * volume;
    }
#ifdef HAVE_MPI
    double tmp(p_error);
    mesh->get_comm()->SumAll(&tmp, &p_error, 1);
    tmp = p_norm;
    mesh->get_comm()->SumAll(&tmp, &p_norm, 1);
#endif

    // calculate flux errors
    Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
    op2->UpdateFlux(solution, flux);
#ifdef HAVE_MPI
    flux.ScatterMasterToGhosted();
#endif

    double flux_norm(0.0), flux_error(0.0);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      const Point& xf = mesh->face_centroid(f);
      const AmanziGeometry::Point& velocity = velocity_exact(xf, 0.0);
      double tmp = velocity * normal;
 
      flux_error += std::pow(tmp - flx[0][f], 2.0);
      flux_norm += std::pow(tmp, 2.0);
    }  
#ifdef HAVE_MPI
    tmp = flux_error;
    mesh->get_comm()->SumAll(&tmp, &flux_error, 1);
    tmp = flux_norm;
    mesh->get_comm()->SumAll(&tmp, &flux_norm, 1);
#endif

    if (MyPID == 0) {
      p_error = pow(p_error / p_norm, 0.5);
      flux_error = pow(flux_error / flux_norm, 0.5);
      printf("scale = %7.4g  Err(p) = %9.6f  Err(flux) = %9.6g  itr=%3d\n", 
          factor, p_error, flux_error, solver->num_itrs()); 
    
      CHECK(p_error< 0.15 && flux_error < 0.15);
    }
  }

  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
  GMV::open_data_file(*mesh, (std::string)"operators.gmv");
  GMV::start_data();
  GMV::write_cell_data(p, 0, "pressure");
  GMV::close_data_file();
}


/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. This is a prototype for future solvers.
* **************************************************************** */
TEST(OPERATOR_NODAL_DIFFUSION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D steady-state elliptic solver, nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
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
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 4, 4, gm);
  RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);
  // RCP<const Mesh> mesh = meshfactory("test/median255x256.exo", gm);

  // create diffusion coefficient
  std::vector<WhetStone::Tensor> K;
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nnodes_owned = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  for (int c = 0; c < ncells_owned; c++) {
    const Point& xc = mesh->cell_centroid(c);
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = Kxx(xc, 0.0);
    Kc(1, 1) = Kyy(xc, 0.0);
    Kc(0, 1) = Kxy(xc, 0.0);
    Kc(1, 0) = Kxy(xc, 0.0);

    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);

  // create boundary data
  Point xv(2);
  std::vector<int> bc_model(nnodes_wghost);
  std::vector<double> bc_value(nnodes_wghost);

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_value[v] = pressure_exact(xv, 0.0);
      bc_model[v] = Operators::OPERATOR_BC_FACE_DIRICHLET;
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(bc_model, bc_value));

  // create diffusion operator 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("node", AmanziMesh::NODE, 1);

  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);

  // create source 
  CompositeVector source(*cvs);
  source.PutScalarMasterAndGhosted(0.0);
  
  Epetra_MultiVector& src = *source.ViewComponent("node", true);
  for (int v = 0; v < nnodes_owned; v++) {
    mesh->node_get_coordinates(v, &xv);
    src[0][v] = source_exact(xv, 0.0);
  }

  // MAIN LOOP
  for (int n = 0; n < 400; n+=110) {
    // double factor = pow(10.0, (double)(n - 50) / 100.0);
    double factor = pow(10.0, (double)(n - 150) / 100.0);

    // create source operator 
    Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(cvs, 0));
    op1->Init();
    op1->UpdateMatrices(source);

    // populate the diffusion operator
    Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operators")
                                        .get<Teuchos::ParameterList>("nodal diffusion");
    Teuchos::RCP<OperatorDiffusion> op2 = Teuchos::rcp(new OperatorDiffusion(*op1, olist, bc));
    int schema_dofs = op2->schema_dofs();
    CHECK(schema_dofs == Operators::OPERATOR_SCHEMA_DOFS_NODE);

    op2->set_factor(factor);  // for developers only
    op2->InitOperator(K, Teuchos::null, Teuchos::null, rho, mu);
    op2->UpdateMatrices(Teuchos::null, Teuchos::null);
    op2->ApplyBCs();
    op2->SymbolicAssembleMatrix(schema_dofs);
    op2->AssembleMatrix(schema_dofs);

    ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
    op2->InitPreconditioner("Hypre AMG", slist);

    // solve the problem
    ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
    solution.PutScalar(0.0);
    AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
       solver = factory.Create("AztecOO CG", lop_list, op2);

    CompositeVector& rhs = *op2->rhs();
    solution.PutScalar(0.0);
    int ierr = solver->ApplyInverse(rhs, solution);

    // calculate errors
#ifdef HAVE_MPI
    solution.ScatterMasterToGhosted();
#endif

    Epetra_MultiVector& sol = *solution.ViewComponent("node", true);
    WhetStone::MFD3D_Diffusion mfd(mesh);
    AmanziGeometry::Point grad(2);

    double p_norm(0.0), p_error(0.0);
    double grad_norm(0.0), grad_error(0.0);
    AmanziMesh::Entity_ID_List nodes;

    for (int c = 0; c < ncells_owned; c++) {
      double volume = mesh->cell_volume(c);

      mesh->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      std::vector<double> cell_solution(nnodes);

      for (int k = 0; k < nnodes; k++) {
        int v = nodes[k];
        cell_solution[k] = sol[0][v];

        mesh->node_get_coordinates(v, &xv);
        double tmp = pressure_exact(xv, 0.0);

        p_error += std::pow(tmp - sol[0][v], 2.0) * volume / nnodes;
        p_norm += std::pow(tmp, 2.0) * volume / nnodes;
      }

      const Point& xc = mesh->cell_centroid(c);
      const AmanziGeometry::Point& grad_exact = gradient_exact(xc, 0.0);
      mfd.RecoverGradient_StiffnessMatrix(c, cell_solution, grad);

      grad_error += L22(grad - grad_exact) * volume;
      grad_norm += L22(grad_exact) * volume;
    }
#ifdef HAVE_MPI
    double tmp(p_error);
    mesh->get_comm()->SumAll(&tmp, &p_error, 1);
    tmp = p_norm;
    mesh->get_comm()->SumAll(&tmp, &p_norm, 1);

    tmp = grad_error;
    mesh->get_comm()->SumAll(&tmp, &grad_error, 1);
    tmp = grad_norm;
    mesh->get_comm()->SumAll(&tmp, &grad_norm, 1);
#endif

    if (MyPID == 0) {
      p_error = pow(p_error / p_norm, 0.5);
      grad_error = pow(grad_error / grad_norm, 0.5);
      double tmp = op2->nfailed_primary() * 100.0 / ncells_owned; 
      printf("scale = %7.4g  Err(p) = %9.6f  Err(grad) = %9.6g  itr=%3d  nfailed=%4.1f\n", 
          factor, p_error, grad_error, solver->num_itrs(), tmp); 

      CHECK(p_error < 0.1 && grad_error < 0.15);
    }
  }
}

