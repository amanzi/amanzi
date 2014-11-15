/*
  This is the operators component of the Amanzi code. 

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

#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
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

double pressure_exact_linear(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];
  return x + 2 * y;
}
Amanzi::AmanziGeometry::Point velocity_exact_linear(const Amanzi::AmanziGeometry::Point& p, double t) { 
  double x = p[0];
  double y = p[1];

  Amanzi::AmanziGeometry::Point v(2);
  v[0] = -5.0;
  v[1] = -3.0;
  return v;
}

int BoundaryFaceGetCell(const Amanzi::AmanziMesh::Mesh& mesh, int f)
{
  Amanzi::AmanziMesh::Entity_ID_List cells;
  mesh.face_get_cells(f, Amanzi::AmanziMesh::USED, &cells);
  return cells[0];
}

/* *****************************************************************
* This test diffusion solver with full tensor and source term.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_NODAL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D steady-state elliptic solver, nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 30, 30, gm);
  // RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = Kxx(xc, 0.0);
    Kc(1, 1) = Kyy(xc, 0.0);
    Kc(0, 1) = Kxy(xc, 0.0);
    Kc(1, 0) = Kxy(xc, 0.0);

    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  Point xv(2);
  std::vector<int> bc_model(nnodes_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nnodes_wghost);

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model[v] = OPERATOR_BC_DIRICHLET;
      bc_value[v] = pressure_exact(xv, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_NODE, bc_model, bc_value));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator nodal");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g);
  const CompositeVectorSpace& cvs = op->DomainMap();
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_NODE;
  op->InitOperator(K, Teuchos::null, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(Operators::OPERATOR_SCHEMA_DOFS_NODE);
  op->AssembleMatrix(schema);

  // create source and add it to the operator
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("node", true);
  src.PutScalar(0.0);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    AmanziMesh::Entity_ID_List nodes;
    mesh->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      src[0][v] += source_exact(xc, 0.0) * volume / nnodes;
    }
  }
  source.GatherGhostedToMaster();

  Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(*op));
  op1->UpdateMatrices(source);

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  Teuchos::RCP<Operator> op2 = Teuchos::rcp(new Operator(*op1));
  op2->InitPreconditioner("Hypre AMG", slist);

  // Test SPD properties of the preconditioner.
  CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);
  a.Random();
  b.Random();
  op2->ApplyInverse(a, ha);
  op2->ApplyInverse(b, hb);

  double ahb, bha, aha, bhb;
  a.Dot(hb, &ahb);
  b.Dot(ha, &bha);
  a.Dot(ha, &aha);
  b.Dot(hb, &bhb);

  if (MyPID == 0) {
    std::cout << "Preconditioner:\n"
              << "  Symmetry test: " << ahb << " = " << bha << std::endl;
    std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
  }
  CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
  CHECK(aha > 0.0);
  CHECK(bhb > 0.0);

  // solve the problem
  ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create("AztecOO CG", lop_list, op2);

  CompositeVector rhs = *op2->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("node");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_node_data(p, 0, "solution");
    GMV::close_data_file();
  }

  CHECK(solver->num_itrs() < 10);
}


/* *****************************************************************
* Exactness test for mixed diffusion solver.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_MIXED) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D steady-state elliptic solver, mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 5, 5, gm);
  RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = 3.0;
    Kc(1, 1) = 1.0;
    Kc(0, 1) = 1.0;
    Kc(1, 0) = 1.0;

    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0), bc_mixed(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    double area = mesh->face_area(f);
    int dir, c = BoundaryFaceGetCell(*mesh, f);
    const Point& normal = mesh->face_normal(f, false, c, &dir);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = velocity_exact_linear(xf, 0.0) * normal / area;  // We assume exterior normal.
    } else if(fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = velocity_exact_linear(xf, 0.0) * normal / area;  // We assume exterior normal.

      double tmp = pressure_exact_linear(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    } else if(fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = pressure_exact_linear(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator mixed");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g);
  const CompositeVectorSpace& cvs = op->DomainMap();
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_FACE;
  op->InitOperator(K, Teuchos::null, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  op->InitPreconditioner("Hypre AMG", slist);

  // Test SPD properties of the preconditioner.
  CompositeVector a(cvs), ha(cvs), b(cvs), hb(cvs);
  a.Random();
  b.Random();
  op->ApplyInverse(a, ha);
  op->ApplyInverse(b, hb);

  double ahb, bha, aha, bhb;
  a.Dot(hb, &ahb);
  b.Dot(ha, &bha);
  a.Dot(ha, &aha);
  b.Dot(hb, &bhb);

  if (MyPID == 0) {
    std::cout << "Preconditioner:\n"
              << "  Symmetry test: " << ahb << " = " << bha << std::endl;
    std::cout << "  Positivity test: " << aha << " " << bhb << std::endl;
  }
  CHECK_CLOSE(ahb, bha, 1e-12 * fabs(ahb));
  CHECK(aha > 0.0);
  CHECK(bhb > 0.0);

  // solve the problem
  ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create("AztecOO CG", lop_list, op);

  CompositeVector rhs = *op->rhs();
  CompositeVector solution(rhs), flux(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;
  }

  // compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double p_norm(0.0), p_error(0.0);
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double tmp = pressure_exact_linear(xc, 0.0);
    double volume = mesh->cell_volume(c);

    // std::cout << c << " " << tmp << " " << p[0][c] << std::endl;
    p_error += std::pow(tmp - p[0][c], 2.0) * volume;
    p_norm += std::pow(tmp, 2.0) * volume;
  }

  // calculate flux error
  Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
  op->UpdateFlux(solution, flux);

  double flux_norm(0.0), flux_error(0.0);
  for (int f = 0; f < nfaces; f++) {
    const Point& normal = mesh->face_normal(f);
    const Point& xf = mesh->face_centroid(f);
    const AmanziGeometry::Point& velocity = velocity_exact_linear(xf, 0.0);
    double tmp = velocity * normal;

    flux_error += std::pow(tmp - flx[0][f], 2.0);
    flux_norm += std::pow(tmp, 2.0);
  }

  if (MyPID == 0) {
    p_error = pow(p_error / p_norm, 0.5);
    flux_error = pow(flux_error / flux_norm, 0.5);
    printf("Err(p) = %9.6f  Err(flux) = %9.6g  itr=%3d\n", p_error, flux_error, solver->num_itrs());

    CHECK(p_error < 1e-12 && flux_error < 1e-12);
    CHECK(solver->num_itrs() < 10);
  }
}


/* *****************************************************************
* Exactness test for mixed diffusion solver.
* NOTE. Mixed boundary condition requires to use mass matrix. We
*       lump it which leads to a small error.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_NODAL_EXACTNESS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D steady-state elliptic solver," 
                            << " exactness test for nodal discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 4, 4, gm);
  RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = 3.0;
    Kc(1, 1) = 1.0;
    Kc(0, 1) = 1.0;
    Kc(1, 0) = 1.0;

    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data.
  Point xv(2);
  std::vector<int> bc_model_v(nnodes_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_v(nnodes_wghost, 0.0);

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if(fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model_v[v] = OPERATOR_BC_DIRICHLET;
      bc_value_v[v] = pressure_exact_linear(xv, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc_v = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_NODE, bc_model_v, bc_value_v));

  std::vector<int> bc_model_f(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_f(nfaces_wghost, 0.0), bc_mixed_f(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6) {
      bc_model_f[f] = OPERATOR_BC_NEUMANN;
      bc_value_f[f] = -(velocity_exact_linear(xf, 0.0))[0];  // We assume exterior normal.
    } else if (fabs(xf[1]) < 1e-6) {
      bc_model_f[f] = OPERATOR_BC_MIXED;
      bc_value_f[f] = -(velocity_exact_linear(xf, 0.0))[1];  // We assume exterior normal.

      double tmp = pressure_exact_linear(xf, 0.0);
      bc_mixed_f[f] = 1.0;
      bc_value_f[f] -= bc_mixed_f[f] * tmp;
    }
  }
  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model_f, bc_value_f, bc_mixed_f));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator nodal");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc_v, op_list, g);
  op->AddBCs(bc_f);
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_NODE;
  op->InitOperator(K, Teuchos::null, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create("AztecOO CG", lop_list, op);

  CompositeVector rhs = *op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;
  }

  // compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("node", false);

  AmanziMesh::Entity_ID_List nodes;
  double p_norm(0.0), p_error(0.0);

  for (int c = 0; c < ncells; c++) {
    mesh->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    double volume = mesh->cell_volume(c);

    for (int n = 0; n < nnodes; n++) {
      int v = nodes[n];
      mesh->node_get_coordinates(v, &xv);
      double tmp = pressure_exact_linear(xv, 0.0);

      // std::cout << v << " " << tmp << " " << p[0][v] << std::endl;
      p_error += std::pow(tmp - p[0][v], 2.0) * volume / nnodes;
      p_norm += std::pow(tmp, 2.0) * volume / nnodes;
    }
  }

  if (MyPID == 0) {
    p_error = pow(p_error / p_norm, 0.5);
    printf("Err(p) = %9.6f itr=%3d\n", p_error, solver->num_itrs());

    CHECK(p_error < 1e-5);
    CHECK(solver->num_itrs() < 10);
  }
}


/* *****************************************************************
* Exactness test for cell-based diffusion solver.
* **************************************************************** */
TEST(OPERATOR_DIFFUSION_CELL_EXACTNESS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: 2D steady-state elliptic solver," 
                            << " exactness test for cell-based discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 5, 8, gm);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    WhetStone::Tensor Kc(2, 2);

    Kc(0, 0) = 3.0;
    Kc(1, 1) = 1.0;
    Kc(0, 1) = 0.0;
    Kc(1, 0) = 0.0;

    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data.
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0), bc_mixed(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    double area = mesh->face_area(f);

    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = 3.0;
    } else if(fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = 2.0;

      double tmp = pressure_exact_linear(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    } else if(fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = pressure_exact_linear(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator cell");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc_f, op_list, g);
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_CELL;
  op->InitOperator(K, Teuchos::null, Teuchos::null, rho, mu);
  op->UpdateMatrices(Teuchos::null, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  // create preconditoner using the base operator class
  ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create("AztecOO CG", lop_list, op);

  CompositeVector rhs = *op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  int ierr = solver->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "pressure solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;
  }

  // compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double p_norm(0.0), p_error(0.0);
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double tmp = pressure_exact_linear(xc, 0.0);
    double volume = mesh->cell_volume(c);

    // std::cout << c << " " << tmp << " " << p[0][c] << std::endl;
    p_error += std::pow(tmp - p[0][c], 2.0) * volume;
    p_norm += std::pow(tmp, 2.0) * volume;
  }

  if (MyPID == 0) {
    p_error = pow(p_error / p_norm, 0.5);
    printf("Err(p) = %9.6f itr=%3d\n", p_error, solver->num_itrs());

    CHECK(p_error < 1e-5);
    CHECK(solver->num_itrs() < 10);
  }
}

