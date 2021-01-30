/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_DIFFUSION_FIXTURE_HH_
#define AMANZI_DIFFUSION_FIXTURE_HH_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

// Amanzi
#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "Verification.hh"
#include "WhetStoneMeshUtils.hh"

#include "AnalyticBase.hh"

using namespace Amanzi;

struct DiffusionFixture {
  DiffusionFixture(Teuchos::RCP<Teuchos::ParameterList> plist_)
    : plist(plist_),
      comm(Amanzi::getDefaultComm()) {};

  void Discretize(const std::string& name, AmanziMesh::Entity_kind kind);
  void DiscretizeWithGravity(const std::string& name, double gravity, AmanziMesh::Entity_kind kind);

  // set up the problem
  void Init(int d, int nx, const std::string& mesh_file);

  // -- general parameters
  void Setup(const std::string& pc_name_, bool symmetric_);

  // -- coefficients
  void SetScalarCoefficient(Operators::PDE_DiffusionFactory& opfactory,
                            AmanziMesh::Entity_kind kind);

  // -- boundary conditions
  void SetBCsDirichlet();
  void SetBCsDirichletNeumann();
  void SetBCsDirichletNeumannRobin();

  void Go(double tol = 1.e-14);
  
  Comm_ptr_type comm;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;

  Teuchos::RCP<AnalyticBase> ana;

  bool symmetric;
  Teuchos::RCP<Operators::BCs> bc;
  Teuchos::RCP<CompositeVector> solution, flux;
  std::string pc_name;
  
  Teuchos::RCP<Operators::PDE_Diffusion> op;
  Teuchos::RCP<Operators::Operator> global_op;
};


/* ******************************************************************
* Initialization
****************************************************************** */
void DiffusionFixture::Init(int d, int nx, const std::string& mesh_file)
{
  if (!plist.get())
    plist = Teuchos::getParametersFromXmlFile("test/operator_diffusion_low_order.xml");

  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(d, plist->sublist("regions"), *comm));
  AmanziMesh::MeshFactory meshfactory(comm, gm);

  if (mesh_file == "structured1d") {
    mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 100, 1);
  } else if (mesh_file == "structured2d") {
    mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, nx, nx);
  } else if (mesh_file == "structured3d") {
    mesh = meshfactory.create(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 4, 5, 6);
  } else {
    mesh = meshfactory.create(mesh_file);
  }
}


/* ******************************************************************
* Create operator
****************************************************************** */
void DiffusionFixture::Discretize(const std::string& name, 
                                  AmanziMesh::Entity_kind scalar_coef)
{
  Operators::PDE_DiffusionFactory opfactory(plist->sublist("PK operator").sublist(name), mesh);

  // populate diffusion coefficient(s)
  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    K->push_back(ana->TensorDiffusivity(xc, 0.0));   
  }

  opfactory.SetVariableTensorCoefficient(K);
  SetScalarCoefficient(opfactory, scalar_coef);

  op = opfactory.Create();

  // boundary condition
  bc = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op->SetBCs(bc, bc);
}


/* ******************************************************************
* Create operator with gravity
****************************************************************** */
void DiffusionFixture::DiscretizeWithGravity(
    const std::string& name, double gravity,
    AmanziMesh::Entity_kind scalar_coef)
{
  Operators::PDE_DiffusionFactory opfactory(plist->sublist("PK operator").sublist(name), mesh);

  // set gravity
  AmanziGeometry::Point g(mesh->space_dimension());
  g[mesh->space_dimension() - 1] = -gravity;
  opfactory.SetConstantGravitationalTerm(g, 1.0);

  // populate diffusion coefficient
  auto K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    K->push_back(ana->TensorDiffusivity(xc, 0.0));   
  }

  opfactory.SetVariableTensorCoefficient(K);
  SetScalarCoefficient(opfactory, scalar_coef);

  op = opfactory.Create();

  bc = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op->SetBCs(bc, bc);
}


/* ******************************************************************
* Set coefficients
****************************************************************** */
void DiffusionFixture::Setup(const std::string& pc_name_, bool symmetric_)
{
  symmetric = symmetric_;
  pc_name = pc_name_;
  global_op = op->global_operator();
  solution = Teuchos::rcp(new CompositeVector(global_op->DomainMap()));
  solution->PutScalar(0.0);
    
  // get and assemble the global operator
  // if (pc_name != "identity") {
  //   global_op->SymbolicAssembleMatrix();
  // }

  // create preconditoner using the base operator class

  if (symmetric) {
    global_op->set_inverse_parameters(pc_name, plist->sublist("preconditioners"),
                                      "AztecOO CG", plist->sublist("solvers"));
    global_op->InitializeInverse();
  } else {
    global_op->set_inverse_parameters(pc_name, plist->sublist("preconditioners"),
                                      "GMRES", plist->sublist("solvers"));
    global_op->InitializeInverse();
  }

  CompositeVectorSpace flux_space;
  flux_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  flux = flux_space.Create();
}


/* ******************************************************************
* Set coefficients
****************************************************************** */
void DiffusionFixture::SetScalarCoefficient(
    Operators::PDE_DiffusionFactory& opfactory,
    AmanziMesh::Entity_kind kind)
{
  if (kind == AmanziMesh::Entity_kind::UNKNOWN) return;
  int nents = mesh->num_entities(kind, AmanziMesh::Parallel_type::ALL);

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)->SetGhosted();
  Teuchos::RCP<CompositeVector> kr;

  if (kind == AmanziMesh::CELL) {
    cvs.SetComponent("cell", AmanziMesh::CELL, 1);
    kr = cvs.Create();
    auto& kr_c = *kr->ViewComponent("cell", true);
    for (int c = 0; c != nents; ++c) {
      kr_c[0][c] = ana->ScalarDiffusivity(mesh->cell_centroid(c), 0.0);
    }

  } else if (kind == AmanziMesh::FACE) {
    cvs.SetComponent("face", AmanziMesh::FACE, 1);
    kr = cvs.Create();
    auto& kr_f = *kr->ViewComponent("face", true);
    for (int f = 0; f != nents; ++f) {
      kr_f[0][f] = ana->ScalarDiffusivity(mesh->face_centroid(f), 0.0);
    }
  }

  opfactory.SetVariableScalarCoefficient(kr);
}


/* ******************************************************************
* Set boundary conditions
****************************************************************** */
void DiffusionFixture::SetBCsDirichlet()
{
  auto& bc_value = bc->bc_value();
  auto& bc_model = bc->bc_model();
    
  if (bc->kind() == AmanziMesh::FACE) {
    const auto& bf_map = mesh->map(AmanziMesh::BOUNDARY_FACE, false);
    const auto& f_map = mesh->map(AmanziMesh::FACE, false);
      
    for (int bf = 0; bf != bf_map.NumMyElements(); ++bf) {
      auto f = f_map.LID(bf_map.GID(bf));
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana->pressure_exact(mesh->face_centroid(f), 0.0);
    }
  } else {
    Exceptions::amanzi_throw("OperatorDiffusion test harness not implemented for this kind.");
  }
}


void DiffusionFixture::SetBCsDirichletNeumann()
{
  auto& bc_value = bc->bc_value();
  auto& bc_model = bc->bc_model();
    
  if (bc->kind() == AmanziMesh::FACE) {
    const auto& bf_map = mesh->map(AmanziMesh::BOUNDARY_FACE, true);
    const auto& f_map = mesh->map(AmanziMesh::FACE, true);
      
    for (int bf = 0; bf != bf_map.NumMyElements(); ++bf) {
      auto f = f_map.LID(bf_map.GID(bf));
      const auto& xf = mesh->face_centroid(f);
      if (xf[0] < 0.0) {
        bool flag;
        auto normal = ana->face_normal_exterior(f, &flag);
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = (ana->velocity_exact(xf, 0.0) * normal) / mesh->face_area(f);
      } else {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana->pressure_exact(xf, 0.0);
      }          
    }
  } else {
    Exceptions::amanzi_throw("OperatorDiffusion test harness not implemented for this kind.");
  }
}


void DiffusionFixture::SetBCsDirichletNeumannRobin()
{
  auto& bc_value = bc->bc_value();
  auto& bc_model = bc->bc_model();
  auto& bc_mixed = bc->bc_mixed();
    
  if (bc->kind() == AmanziMesh::FACE) {
    const auto& bf_map = mesh->map(AmanziMesh::BOUNDARY_FACE, true);
    const auto& f_map = mesh->map(AmanziMesh::FACE, true);
      
    for (int bf = 0; bf != bf_map.NumMyElements(); ++bf) {
      auto f = f_map.LID(bf_map.GID(bf));
      const auto& xf = mesh->face_centroid(f);

      bool flag;
      auto normal = ana->face_normal_exterior(f, &flag);
      double area = mesh->face_area(f);

      if (xf[0] < 0.0) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = (ana->velocity_exact(xf, 0.0) * normal) / mesh->face_area(f);
      } else if (xf[1] < 0.0) {
        bc_model[f] = Operators::OPERATOR_BC_MIXED;
        bc_value[f] = ana->velocity_exact(xf, 0.0) * normal / area;

        double tmp = ana->pressure_exact(xf, 0.0);
        bc_mixed[f] = 1.0;
        bc_value[f] -= bc_mixed[f] * tmp;
      } else {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana->pressure_exact(xf, 0.0);
      }          
    }
  } else {
    Exceptions::amanzi_throw("OperatorDiffusion test harness not implemented for this kind.");
  }
}


/* ******************************************************************
* Run the solvers
****************************************************************** */
void DiffusionFixture::Go(double tol)
{
  global_op->Init();
  op->UpdateMatrices(Teuchos::null, solution.ptr());

  int ncells = mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  CompositeVector& rhs = *global_op->rhs();

  auto& rhs_c = *rhs.ViewComponent("cell", false);
  for (int c = 0; c != ncells; ++c) {
    const auto& xc = mesh->cell_centroid(c);
    rhs_c[0][c] += ana->source_exact(xc, 0.0) * mesh->cell_volume(c);
  }
    
  op->ApplyBCs(true, true, true);
  global_op->ComputeInverse();
  global_op->ApplyInverse(*global_op->rhs(), *solution);

  auto MyPID = comm->MyPID();
  if (MyPID == 0) {
    std::cout << "pressure solve: ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() << std::endl;
  }
      
  if (tol > 0.0) {
    // compute pressure error
    Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
    double pnorm(0.0), pl2_err(0.0), pinf_err(0.0);
    ana->ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);

    // calculate flux error
    double unorm, ul2_err, uinf_err;

    op->UpdateMatrices(Teuchos::null, solution.ptr());
    op->UpdateFlux(solution.ptr(), flux.ptr());
    {
      auto fv = flux->ViewComponent("face", false);
    }

    Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
    ana->ComputeFaceError(flx, 0.0, unorm, ul2_err, uinf_err);

    if (MyPID == 0) {
      pl2_err /= pnorm; 
      ul2_err /= unorm;
      printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n",
             pl2_err, pinf_err, ul2_err, uinf_err, global_op->num_itrs());
       
      CHECK(pl2_err < tol);
      CHECK(ul2_err < 10*tol);
      //if (pc_name != "identity" && pc_name != "diagonal") CHECK(solver->num_itrs() < 10);
    }

    if (symmetric) {
      VerificationCV ver(global_op);
      ver.CheckMatrixSPD();
    }
  }
}

#endif
