/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

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
  ParameterList parameter_list;
  string xmlFileName = "test/flow_darcy_source.xml";
  
  // DEPRECATED  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 200, 200, gm);
  RCP<Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // create and populate fake flow state
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(mesh));
  FS->Initialize();
  FS->set_permeability_2D(2.0, 1.0, "Material 1");
  FS->set_porosity(1.0);
  FS->set_fluid_viscosity(1.0);
  FS->set_fluid_density(1.0);
  FS->set_gravity(0.0);

  // create Darcy process kernel
  Darcy_PK* DPK = new Darcy_PK(parameter_list, FS);
  DPK->InitPK();
  DPK->InitSteadyState(0.0, 1e-8);

  // reset problem coefficients
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  std::vector<WhetStone::Tensor>& K = DPK->get_K();

  for (int c = 0; c < K.size(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    K[c](0, 0) = Kxx(xc, 0.0);
    K[c](1, 1) = Kyy(xc, 0.0);
    K[c](0, 1) = Kxy(xc, 0.0);
    K[c](1, 0) = Kxy(xc, 0.0);
  }

  // update Dirichlet boundary data 
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
    DPK->matrix()->CreateMFDmassMatrices_ScaledStability(
        AmanziFlow::FLOW_MFD3D_POLYHEDRA, factor, K);

    // create a problem
    DPK->AssembleMatrixMFD();

    // update right-hand side
    Teuchos::RCP<Epetra_Vector> rhs = DPK->matrix()->rhs();
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      double volume = mesh->cell_volume(c);
      (*rhs)[c] += source_exact(xc, 0.0) * volume;
    }

    // steady-state solution
    Epetra_Vector& solution = DPK->ref_solution();
    DPK->SolveFullySaturatedProblem(0.0, *rhs, solution);
    DPK->CommitState(FS);

    // calculate errors
    Epetra_Vector& p = FS->ref_pressure();
    Epetra_Vector& flux = FS->ref_darcy_flux();

    double p_norm(0.0), p_error(0.0);
    for (int c = 0; c < ncells; c++) {
      const Point& xc = mesh->cell_centroid(c);
      double tmp = pressure_exact(xc, 0.0);
      double volume = mesh->cell_volume(c);

      p_error += std::pow(tmp - p[c], 2.0) * volume;
      p_norm += std::pow(tmp, 2.0) * volume;
    }

    double flux_norm(0.0), flux_error(0.0);
    for (int f = 0; f < nfaces; f++) {
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      const Point& xf = mesh->face_centroid(f);
      const AmanziGeometry::Point& velocity = velocity_exact(xf, 0.0);
      double tmp = velocity * normal;
 
      flux_error += std::pow(tmp - flux[f], 2.0);
      flux_norm += std::pow(tmp, 2.0);
    }  

    p_error = pow(p_error / p_norm, 0.5);
    flux_error = pow(flux_error / flux_norm, 0.5);
    printf("scale = %10.6g  Err(p) = %10.6f  Err(flux) = %10.6g\n", factor, p_error, flux_error); 
    
    CHECK(p_error< 0.15 && flux_error < 0.15);
  }

  Epetra_Vector& p = FS->ref_pressure();
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
  ParameterList parameter_list;
  string xmlFileName = "test/flow_darcy_source.xml";
  
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 80, 80, gm);
  RCP<Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // create and populate fake flow state
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(mesh));
  FS->Initialize();
  FS->set_permeability_2D(2.0, 1.0, "Material 1");
  FS->set_porosity(1.0);
  FS->set_fluid_viscosity(1.0);
  FS->set_fluid_density(1.0);
  FS->set_gravity(0.0);

  // create diffusion coefficient
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
  const Epetra_Map& map = mesh->node_map(false);
  Stiffness_MFD matrix(FS, map);
  Stiffness_MFD preconditioner(FS, map);

  matrix.SymbolicAssembleGlobalMatrices();
  preconditioner.SymbolicAssembleGlobalMatrices();

  for (int n = 0; n < 400; n+=100) {
    // double factor = pow(10.0, (double)(n - 50) / 100.0);
    double factor = pow(10.0, (double)(n - 150) / 100.0);

    // populate matrix
    matrix.CreateMFDstiffnessMatrices(K, factor);
    matrix.CreateMFDrhsVectors();
    matrix.ApplyBoundaryConditions(bc_model, bc_values);
    matrix.AssembleGlobalMatrices();

    // populate preconditioner
    preconditioner.CreateMFDstiffnessMatrices(K, factor);
    preconditioner.CreateMFDrhsVectors();
    preconditioner.ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner.AssembleGlobalMatrices();
    preconditioner.UpdatePreconditioner();
  
    // update right-hand side
    AmanziMesh::Entity_ID_List nodes;
    Teuchos::RCP<Epetra_Vector> rhs = matrix.rhs();

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
    AztecOO* solver = new AztecOO;
    solver->SetUserOperator(&matrix);
    solver->SetPrecOperator(&preconditioner);

    solver->SetAztecOption(AZ_solver, AZ_cg);
    solver->SetAztecOption(AZ_output, 0);
    solver->SetAztecOption(AZ_conv, AZ_rhs);

    Epetra_Vector solution(map);
    solver->SetRHS(&*rhs);
    solver->SetLHS(&solution); 

    solver->Iterate(100, 1e-16);
    delete solver;

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
