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
#include "UpwindSecondOrder.hh"
#include "UpwindStandard.hh"

#include "Analytic01.hh"
#include "Analytic02.hh"
#include "Analytic03.hh"

namespace Amanzi{

// This class wraps scalar diffusion coefficient Analytic03.
class HeatConduction {
 public:
  HeatConduction(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh), ana_(mesh) { 
    int dim = mesh_->space_dimension();
    cvs_.SetMesh(mesh_);
    cvs_.SetGhosted(true);
    cvs_.SetComponent("cell", AmanziMesh::CELL, 1);
    cvs_.SetOwned(false);
    cvs_.AddComponent("face", AmanziMesh::FACE, 1);
    cvs_.AddComponent("grad", AmanziMesh::CELL, dim);

    values_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
    derivatives_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
  }
  ~HeatConduction() {};

  // main members
  void UpdateValues(const CompositeVector& u) { 
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true); 
    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      const WhetStone::Tensor& Kc = ana_.Tensor(xc, 0.0);
      vcell[0][c] = Kc(0, 0);
    }

    // add gradient component
    int dim = mesh_->space_dimension();
    Epetra_MultiVector& vgrad = *values_->ViewComponent("grad", true); 

    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      AmanziGeometry::Point grad = ana_.ScalarTensorGradient(xc, 0.0);
      for (int i = 0; i < dim; i++) vgrad[i][c] = grad[i];
    }

    derivatives_->PutScalar(1.0);
  }

  void UpdateValuesPostUpwind() { 
    if (!values_->HasComponent("twin")) {
      cvs_.AddComponent("twin", AmanziMesh::FACE, 1);
      Teuchos::RCP<CompositeVector> tmp = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));

      *tmp->ViewComponent("cell") = *values_->ViewComponent("cell"); 
      *tmp->ViewComponent("face") = *values_->ViewComponent("face"); 
      *tmp->ViewComponent("grad") = *values_->ViewComponent("grad"); 
      values_ = tmp;
    }

    AmanziMesh::Entity_ID_List cells;
    Epetra_MultiVector& vcell = *values_->ViewComponent("cell", true); 
    Epetra_MultiVector& vface = *values_->ViewComponent("face", true); 
    Epetra_MultiVector& vtwin = *values_->ViewComponent("twin", true); 

    vtwin = vface;
    int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

    for (int f = 0; f < nfaces; f++) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();
      
      if (ncells == 2) {
        double v1 = vcell[0][cells[0]];
        double v2 = vcell[0][cells[1]];
        if (fabs(v1 - v2) > 2 * std::min(fabs(v1), fabs(v2))) {
          vface[0][f] = v1;
          vtwin[0][f] = v2;
        }  
      } 
    }
  }

  double Conduction(int c, double T) const {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana_.Tensor(xc, 0.0);
    return Kc(0, 0);
  }

  Teuchos::RCP<CompositeVector> values() { return values_; }
  Teuchos::RCP<CompositeVector> derivatives() { return derivatives_; }
   
 private:
  CompositeVectorSpace cvs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<CompositeVector> values_, derivatives_;
  mutable Analytic03 ana_;
};

typedef double(HeatConduction::*ModelUpwindFn)(int c, double T) const; 
}  // namespace Amanzi


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

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, nodal discretization" << std::endl;

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

  Analytic01 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  Point xv(2);
  std::vector<int> bc_model(nnodes_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nnodes_wghost);
  std::vector<double> bc_mixed;

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if (fabs(xv[0]) < 1e-6 || fabs(xv[0] - 1.0) < 1e-6 ||
        fabs(xv[1]) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model[v] = OPERATOR_BC_DIRICHLET;
      bc_value[v] = ana.pressure_exact(xv, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_NODE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator nodal");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g, 0);
  const CompositeVectorSpace& cvs = op->DomainMap();
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_NODE;
  op->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
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
      src[0][v] += ana.source_exact(xc, 0.0) * volume / nnodes;
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
* Tests DivK diffusion solver with full tensor and source term.
* The model for kf is volime-weighted arithmetic average.
***************************************************************** */
TEST(OPERATOR_DIFFUSION_DIVK_AVERAGE_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, divK discretization, average" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator divk");

  // create an SIMPLE mesh framework
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, NULL);
  std::string file = op_list.get<std::string>("file name", "test/random20.exo");
  Teuchos::RCP<const Mesh> mesh = meshfactory(file, NULL);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Analytic03 ana(mesh);

  const WhetStone::Tensor Kc(2, 1);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells; c++) K.push_back(Kc);

  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0), bc_mixed(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g, 0);
  const CompositeVectorSpace& cvs = op->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  Point velocity(-1.0, 0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& normal = mesh->face_normal(f);
    flx[0][f] = velocity * normal;
  }
  
  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh));

  // Create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
  UpwindStandard<HeatConduction> upwind(mesh, knc);
  upwind.Init(ulist);

  knc->UpdateValues(*flux);  // argument is not used
  ModelUpwindFn func = &HeatConduction::Conduction;
  upwind.Compute(*flux, bc_model, bc_value, *knc->values(), *knc->values(), func);

  // knc->UpdateValuesPostUpwind();

  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_FACE;
  op->Setup(K, knc->values(), knc->derivatives(), rho, mu);
  op->UpdateMatrices(flux, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // dirty trick: add source to operator
  Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(*op));
  op1->UpdateMatrices(source);

  // create preconditoner using the base operator class
  Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  Teuchos::ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
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
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error
  op->UpdateFlux(solution, *flux);
  double unorm, ul2_err, uinf_err;
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%12.8g  Inf(p)=%12.8g  L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, solver->num_itrs());

    CHECK(pl2_err < 0.03 && ul2_err < 0.1);
    CHECK(solver->num_itrs() < 10);
  }
}


/* *****************************************************************
* Tests DivK diffusion solver with full tensor and source term.
* The model for kf is volime-weighted arithmetic average.
***************************************************************** */
TEST(OPERATOR_DIFFUSION_DIVK_AVERAGE_3D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 3D elliptic solver, divK discretization, average" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10, NULL);
  // Teuchos::RCP<const Mesh> mesh = meshfactory("test/mesh.exo", NULL);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Analytic03 ana(mesh);

  const WhetStone::Tensor Kc(3, 1);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells; c++) K.push_back(Kc);

  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, 0.0, -1.0);

  // create boundary data
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0), bc_mixed(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  Teuchos::ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator divk");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g, 0);
  const CompositeVectorSpace& cvs = op->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  Point velocity(0.0, 0.0, 0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& normal = mesh->face_normal(f);
    flx[0][f] = velocity * normal;
  }
  
  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh));

  // Create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
  UpwindStandard<HeatConduction> upwind(mesh, knc);
  upwind.Init(ulist);

  knc->UpdateValues(*flux);  // argument is not used
  ModelUpwindFn func = &HeatConduction::Conduction;
  upwind.Compute(*flux, bc_model, bc_value, *knc->values(), *knc->values(), func);

  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_FACE;
  op->Setup(K, knc->values(), knc->derivatives(), rho, mu);
  op->UpdateMatrices(flux, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // dirty trick: add source to operator
  Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(*op));
  op1->UpdateMatrices(source);

  // create preconditoner using the base operator class
  Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  Teuchos::ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
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
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error
  op->UpdateFlux(solution, *flux);
  double unorm, ul2_err, uinf_err;
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, solver->num_itrs());

    CHECK(pl2_err < 0.03 && ul2_err < 0.1);
    CHECK(solver->num_itrs() < 10);
  }
}


/* *****************************************************************
* Tests DivK diffusion solver with full tensor and source term.
* The model for kf is second-order upwind with arithmetic average.
***************************************************************** */
TEST(OPERATOR_DIFFUSION_SECOND_ORDER) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, divK discretization, 2nd-order" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator")
                                        .sublist("diffusion operator second-order");

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  // Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 5, 5, gm);
  std::string file = op_list.get<std::string>("file name", "test/random20.exo");
  Teuchos::RCP<const Mesh> mesh = meshfactory(file, NULL);

  /* modify diffusion coefficient */
  std::vector<WhetStone::Tensor> K;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Analytic03 ana(mesh);

  const WhetStone::Tensor Kc(2, 1);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells; c++) {
    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0), bc_mixed(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g, 0);
  const CompositeVectorSpace& cvs = op->DomainMap();

  // create and initialize state variables.
  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(cvs));
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);

  Point velocity(-1.0, 0.0);
  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& normal = mesh->face_normal(f);
    flx[0][f] = velocity * normal;
  }
  
  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh));

  // Create upwind model
  Teuchos::ParameterList& ulist = plist.sublist("PK operator").sublist("upwind second-order");
  UpwindSecondOrder<HeatConduction> upwind(mesh, knc);
  upwind.Init(ulist);

  knc->UpdateValues(*flux);  // argument is not used
  ModelUpwindFn func = &HeatConduction::Conduction;
  upwind.Compute(*flux, bc_model, bc_value, *knc->values(), *knc->values(), func);

  knc->UpdateValuesPostUpwind();

  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_FACE;
  op->Setup(K, knc->values(), knc->derivatives(), rho, mu);
  op->UpdateMatrices(flux, Teuchos::null);
  op->ApplyBCs();
  op->SymbolicAssembleMatrix(schema);
  op->AssembleMatrix(schema);

  // create source 
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0);
  }

  // dirty trick: add source to operator
  Teuchos::RCP<OperatorSource> op1 = Teuchos::rcp(new OperatorSource(*op));
  op1->UpdateMatrices(source);

  // create preconditoner using the base operator class
  Teuchos::ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  op->InitPreconditioner("Hypre AMG", slist);

  // solve the problem
  Teuchos::ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
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
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error
  op->UpdateFlux(solution, *flux);
  double unorm, ul2_err, uinf_err;
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%12.8g  Inf(p)=%12.8g  L2(u)=%12.8g  Inf(u)=%12.8g  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, solver->num_itrs());

    CHECK(pl2_err < 0.03 && ul2_err < 0.1);
    CHECK(solver->num_itrs() < 10);
  }
}


/* *****************************************************************
* Exactness test for mixed diffusion solver.
***************************************************************** */
TEST(OPERATOR_DIFFUSION_MIXED) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
                            << " test for mixed discretization" << std::endl;

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

  Analytic02 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
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
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;  // We assume exterior normal.
    } else if(fabs(xf[1]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_MIXED;
      bc_value[f] = ana.velocity_exact(xf, 0.0) * normal / area;  // We assume exterior normal.

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    } else if(fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator mixed");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc, op_list, g, 0);
  const CompositeVectorSpace& cvs = op->DomainMap();
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_FACE;
  op->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
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
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  // calculate flux error
  Epetra_MultiVector& flx = *flux.ViewComponent("face", true);
  double unorm, ul2_err, uinf_err;

  op->UpdateFlux(solution, flux);
  ana.ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    ul2_err /= unorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f  itr=%3d\n",
        pl2_err, pinf_err, ul2_err, uinf_err, solver->num_itrs());

    CHECK(pl2_err < 1e-12 && ul2_err < 1e-12);
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
  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
                            << " test for nodal discretization" << std::endl;

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

  Analytic02 ana(mesh);

  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    const WhetStone::Tensor& Kc = ana.Tensor(xc, 0.0);
    K.push_back(Kc);
  }
  double rho(1.0), mu(1.0);
  AmanziGeometry::Point g(0.0, -1.0);

  // create boundary data (no mixed bc)
  Point xv(2);
  std::vector<int> bc_model_v(nnodes_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_v(nnodes_wghost, 0.0);
  std::vector<double> bc_mixed_v;

  for (int v = 0; v < nnodes_wghost; v++) {
    mesh->node_get_coordinates(v, &xv);
    if(fabs(xv[0] - 1.0) < 1e-6 || fabs(xv[1] - 1.0) < 1e-6) {
      bc_model_v[v] = OPERATOR_BC_DIRICHLET;
      bc_value_v[v] = ana.pressure_exact(xv, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc_v = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_NODE, bc_model_v, bc_value_v, bc_mixed_v));

  std::vector<int> bc_model_f(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_f(nfaces_wghost, 0.0), bc_mixed_f(nfaces_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6) {
      bc_model_f[f] = OPERATOR_BC_NEUMANN;
      bc_value_f[f] = -(ana.velocity_exact(xf, 0.0))[0];  // We assume exterior normal.
    } else if (fabs(xf[1]) < 1e-6) {
      bc_model_f[f] = OPERATOR_BC_MIXED;
      bc_value_f[f] = -(ana.velocity_exact(xf, 0.0))[1];  // We assume exterior normal.

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed_f[f] = 1.0;
      bc_value_f[f] -= bc_mixed_f[f] * tmp;
    }
  }
  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model_f, bc_value_f, bc_mixed_f));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator nodal");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc_v, op_list, g, 0);
  op->AddBCs(bc_f);
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_NODE;
  op->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
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

  double pnorm, pl2_err, pinf_err, hnorm, ph1_err;
  ana.ComputeNodeError(p, 0.0, pnorm, pl2_err, pinf_err, hnorm, ph1_err);

  if (MyPID == 0) {
    pl2_err /= pnorm;
    ph1_err /= hnorm;
    printf("L2(p)=%9.6f  H1(p)=%9.6f  itr=%3d\n", pl2_err, ph1_err, solver->num_itrs());

    CHECK(pl2_err < 1e-5 && ph1_err < 2e-5);
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

  if (MyPID == 0) std::cout << "\nTest: 2D elliptic solver, exactness" 
                            << " test for cell-based discretization" << std::endl;

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

  Analytic02 ana(mesh);

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

      double tmp = ana.pressure_exact(xf, 0.0);
      bc_mixed[f] = 1.0;
      bc_value[f] -= bc_mixed[f] * tmp;
    } else if(fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana.pressure_exact(xf, 0.0);
    }
  }
  Teuchos::RCP<BCs> bc_f = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // create diffusion operator 
  ParameterList op_list = plist.get<Teuchos::ParameterList>("PK operator").sublist("diffusion operator cell");
  OperatorDiffusionFactory opfactory;
  Teuchos::RCP<OperatorDiffusion> op = opfactory.Create(mesh, bc_f, op_list, g, 0);
  
  // populate the diffusion operator
  int schema = Operators::OPERATOR_SCHEMA_DOFS_CELL;
  op->Setup(K, Teuchos::null, Teuchos::null, rho, mu);
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
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

  if (MyPID == 0) {
    pl2_err /= pnorm; 
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, solver->num_itrs());

    CHECK(pl2_err < 1e-5);
    CHECK(solver->num_itrs() < 10);
  }
}

