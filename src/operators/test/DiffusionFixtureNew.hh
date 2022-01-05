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
// #include "WhetStoneMeshUtils.hh"

#include "AnalyticBase.hh"

using namespace Amanzi;

struct DiffusionFixture {
  DiffusionFixture(Teuchos::RCP<Teuchos::ParameterList> plist_)
    : plist(plist_),
      comm(Amanzi::getDefaultComm()) {};

  template<class PDE_Diffusion_type>
  void Discretize(const std::string& name, AmanziMesh::Entity_kind kind);

  // set up the problem
  void Init(int d, int nx, const std::string& mesh_file);

  // -- general parameters
  void Setup(const std::string& prec_solver_, bool symmetric_);

  // -- coefficients
  void SetScalarCoefficient(Operators::PDE_DiffusionFactory& opfactory,
                            AmanziMesh::Entity_kind kind);

  // -- boundary conditions
  void SetBCsDirichlet();
  void SetBCsDirichletNeumann() {};
  void SetBCsDirichletNeumannRobin() {};

  // main loop:
  //   tolerance of a Krylov solver
  //   if initial_guess=true, solution from last iteration is used as the initial guess
  void Go(double tol = 1.0e-12, bool initial_guess = true);
  void MatVec(int nloops);

  // access
  Comm_ptr_type get_comm() const { return comm; }
  
 public:
  Comm_ptr_type comm;
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;

  Teuchos::RCP<AnalyticBase> ana;

  bool symmetric;
  Teuchos::RCP<Operators::BCs> bc;
  Teuchos::RCP<CompositeVector> solution, flux;
  std::string prec_solver;
  
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
    mesh = meshfactory.create(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, nx, nx, nx);
  } else {
    mesh = meshfactory.create(mesh_file);
  }
}


/* ******************************************************************
* Create operator
****************************************************************** */
template<class PDE_Diffusion_type>
void DiffusionFixture::Discretize(const std::string& name, 
                                  AmanziMesh::Entity_kind scalar_coef)
{
  op = Teuchos::rcp(new PDE_Diffusion_type(plist->sublist("PK operator").sublist(name), mesh));
  op->Init();
  // modify diffusion coefficient
  CompositeVectorSpace K_map;
  K_map.SetMesh(mesh);
  K_map.AddComponent("cell", AmanziMesh::CELL, 1);
  auto K = Teuchos::rcp(new TensorVector(K_map));
  K->Init(
      K->size(), 
      //size function: size of element c 
      [&](int c) -> const Amanzi::WhetStone::Tensor<Kokkos::HostSpace>& {
       const AmanziGeometry::Point& xc = mesh->cell_centroid_host(c); 
       return ana->TensorDiffusivity(xc,0.0); 
      }
  );
  op->SetTensorCoefficient(K);

  if (scalar_coef == AmanziMesh::Entity_kind::FACE) {
    int nents = mesh->num_entities(scalar_coef, AmanziMesh::Parallel_type::ALL);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh)->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
    Teuchos::RCP<CompositeVector> kr = cvs.Create();
    auto vec = kr->ViewComponent<MirrorHost>("face", true);
    for (int f = 0; f != nents; ++f) {
      vec(f, 0) = ana->ScalarDiffusivity(mesh->face_centroid_host(f), 0.0);
    }
    op->SetScalarCoefficient(kr, Teuchos::null);
  }

  // boundary condition
  bc = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op->SetBCs(bc,bc);
}


/* ******************************************************************
* Set coefficients
****************************************************************** */
void DiffusionFixture::Setup(const std::string& prec_solver_, bool symmetric_)
{
  symmetric = symmetric_;
  prec_solver = prec_solver_;
  global_op = op->global_operator();
  solution = Teuchos::rcp(new CompositeVector(global_op->getDomainMap()));
  solution->putScalar(0.0);
    
  // create preconditoner using the base operator class
  if (prec_solver.substr(0, 6) == "Amesos") {
    global_op->set_inverse_parameters(prec_solver, plist->sublist("solvers"));
  } else {
    if (symmetric) {
      global_op->set_inverse_parameters(prec_solver, plist->sublist("preconditioners"),
                                        "AztecOO CG", plist->sublist("solvers"));
    } else {
      global_op->set_inverse_parameters(prec_solver, plist->sublist("preconditioners"),
                                      "GMRES", plist->sublist("solvers"));
    } 
  }
  global_op->initializeInverse();

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
    auto vec = kr->ViewComponent<MirrorHost>("cell", true);
    for (int c = 0; c != nents; ++c) {
      vec(c, 0) = ana->ScalarDiffusivity(mesh->cell_centroid_host(c), 0.0);
    }

  } else if (kind == AmanziMesh::FACE) {
    cvs.SetComponent("face", AmanziMesh::FACE, 1);
    kr = cvs.Create();
    auto vec = kr->ViewComponent<MirrorHost>("face", true);
    for (int f = 0; f != nents; ++f) {
      vec(f, 0) = ana->ScalarDiffusivity(mesh->face_centroid_host(f), 0.0);
    }
  }
  op->SetScalarCoefficient(kr, Teuchos::null);
}


/* ******************************************************************
* Set boundary conditions
****************************************************************** */
void DiffusionFixture::SetBCsDirichlet()
{
  auto bc_value = bc->bc_value();
  auto bc_model = bc->bc_model();
    
  if (bc->kind() == AmanziMesh::FACE) {
    const auto& bf_map = *mesh->map(AmanziMesh::BOUNDARY_FACE, false);
    const auto& f_map = *mesh->map(AmanziMesh::FACE, false);
      
    for (int bf = 0; bf != bf_map.getNodeNumElements(); ++bf) {
      auto f = f_map.getLocalElement(bf_map.getGlobalElement(bf));
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = ana->pressure_exact(mesh->face_centroid_host(f), 0.0);
    }
  } else {
    Exceptions::amanzi_throw("OperatorDiffusion test harness not implemented for this kind.");
  }
}


/* ******************************************************************
* Run the solvers
****************************************************************** */
void DiffusionFixture::Go(double tol, bool initial_guess)
{
  global_op->Zero();
  op->UpdateMatrices(Teuchos::null, solution.ptr());
  int ncells = mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  CompositeVector& rhs = *global_op->rhs();

  // The view on host needs to be scoped to be released 
  // This would prevent the vector to be copied later (in apply inverse)
  {
    auto rhs_c = rhs.ViewComponent<MirrorHost>("cell", false);
    for (int c = 0; c != ncells; ++c) {
      const auto& xc = mesh->cell_centroid_host(c);
      rhs_c(c, 0) += ana->source_exact(xc, 0.0) * mesh->cell_volume_host(c);
      std::cout<<rhs_c(c,0)<<"("<<ana->source_exact(xc, 0.0)<<"*"<<mesh->cell_volume_host(c)<<") - ";
    }
    std::cout<<std::endl;
  }

  

  op->ApplyBCs(true, true, true);
  global_op->computeInverse();
  if (!initial_guess) solution->putScalar(0.0);
  auto start = std::chrono::high_resolution_clock::now();
  global_op->applyInverse(rhs, *solution);
  auto stop = std::chrono::high_resolution_clock::now();
  if (tol > 0.0) {
    // compute pressure error
    double pnorm(0.0), pl2_err(0.0), pinf_err(0.0);
    ComputeCellError(*ana, mesh, *solution, 0.0, pnorm, pl2_err, pinf_err);

    // calculate flux error
    op->UpdateMatrices(Teuchos::null, solution.ptr());
    op->UpdateFlux(solution.ptr(), flux.ptr());

    double unorm, ul2_err, uinf_err;
    ComputeFaceError(*ana, mesh, *flux, 0.0, unorm, ul2_err, uinf_err);

    auto MyPID = comm->getRank();
    if (MyPID == 0) {
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      printf("iterative loop time = %g sec\n", double(duration.count() * 1e-6));

      pl2_err /= pnorm; 
      ul2_err /= unorm;
      printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f\ndofs=    itr=%d  ||r||=%10.5f\n",
             pl2_err, pinf_err, ul2_err, uinf_err,
             global_op->num_itrs(), global_op->residual());
       
      CHECK(pl2_err < 10*tol);
      CHECK(ul2_err < 100*tol);
    }
    // if (symmetric) {
    //   VerificationCV ver(global_op);
    //   ver.CheckMatrixSPD();
    // }
  }
}


/* ******************************************************************
* Run the mat-vec
****************************************************************** */
void DiffusionFixture::MatVec(int nloops)
{
  global_op->Zero();
  op->UpdateMatrices(Teuchos::null, solution.ptr());
  int ncells = mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

  op->ApplyBCs(true, true, true);
  global_op->computeInverse();

  auto sol(*solution);

  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < nloops; ++i) {
    global_op->apply(sol, *solution);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  printf("iterative loop time = %g sec\n", double(duration.count() * 1e-6));
}

#endif
