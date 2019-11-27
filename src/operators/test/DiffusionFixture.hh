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
#include "Teuchos_ParameterXMLFileReader.hpp"

// Amanzi
#include "MeshFactory.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "Tensor.hh"


#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"

using namespace Amanzi;

struct DiffusionFixture {
  DiffusionFixture(const Teuchos::RCP<AnalyticBase>& ana_,
                   const std::string& mesh_file)
      : comm(Amanzi::getDefaultComm()),
        ana(ana_)
  {
    std::string xmlFileName = "test/operator_diffusion_discretizations.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    plist = xmlreader.getParameters();

    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(ana->dimension(),
            plist.sublist("regions"), *comm));
    AmanziMesh::MeshFactory meshfactory(comm, gm);
    if (mesh_file == "Generate2D") {
      mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 10, 10);
    } else if (mesh_file == "Generate3D") {
      mesh = meshfactory.create(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 4, 5, 6);
    } else {
      mesh = meshfactory.create(mesh_file);
    }
    
  }

  template<class PDE_Diffusion_type, AmanziMesh::Entity_kind Boundary_kind>
  void discretize(const std::string& name) {
    // construct the operator
    op = Teuchos::rcp(new PDE_Diffusion_type(plist.sublist("PK operator").sublist(name), mesh));
    op->Init();

    // modify diffusion coefficient
    Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      const WhetStone::Tensor Kc = ana->TensorDiffusivity(xc, 0.0);
      K->push_back(Kc);
    }
    op->SetTensorCoefficient(K);

    // boundary condition
    bc = Teuchos::rcp(new Operators::BCs(mesh, Boundary_kind, WhetStone::DOF_Type::SCALAR));
    op->SetBCs(bc,bc);
  }

  template<class PDE_Diffusion_type, AmanziMesh::Entity_kind Boundary_kind>
  void discretizeWithGravity(const std::string& name, double gravity) {
    AmanziGeometry::Point g(mesh->space_dimension());
    g[mesh->space_dimension()-1] = -gravity;
    Teuchos::ParameterList op_list = plist.sublist("PK operator").sublist(name);
    op_list.set("gravity", true);
    op = Teuchos::rcp(new PDE_Diffusion_type(op_list, mesh));
    op->SetDensity(1.0);
    op->SetGravity(g);
    op->Init();

    // modify diffusion coefficient
    Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c < ncells; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      const WhetStone::Tensor& Kc = ana->TensorDiffusivity(xc, 0.0);
      K->push_back(Kc);
    }
    op->SetTensorCoefficient(K);

    // boundary condition
    bc = Teuchos::rcp(new Operators::BCs(mesh, Boundary_kind, WhetStone::DOF_Type::SCALAR));
    op->SetBCs(bc,bc);
  }

  void scalarCoefficient(AmanziMesh::Entity_kind kind) {
    CompositeVectorSpace space;
    space.SetMesh(mesh)->SetGhosted();
    Teuchos::RCP<CompositeVector> kr;
    if (kind == AmanziMesh::CELL) {
      space.SetComponent("cell", AmanziMesh::CELL, 1);
      kr = space.Create();
      auto& vec = *kr->ViewComponent("cell", true);
      for (int i=0; i!=mesh->num_entities(kind, AmanziMesh::Parallel_type::ALL); ++i) {
        vec[0][i] = ana->ScalarDiffusivity(mesh->cell_centroid(i), 0.0);
      }

    } else if (kind == AmanziMesh::FACE) {
      space.SetComponent("face", AmanziMesh::FACE, 1);
      kr = space.Create();
      auto& vec = *kr->ViewComponent("face", true);
      for (int i=0; i!=mesh->num_entities(kind, AmanziMesh::Parallel_type::ALL); ++i) {
        vec[0][i] = ana->ScalarDiffusivity(mesh->face_centroid(i), 0.0);
      }
    }
    op->SetScalarCoefficient(kr, Teuchos::null);
  }

  void setBCsDirichlet() {
    auto& bc_value = bc->bc_value();
    auto& bc_model = bc->bc_model();
    
    if (bc->kind() == AmanziMesh::FACE) {
      auto bf_map = mesh->map(AmanziMesh::BOUNDARY_FACE, false);
      auto f_map = mesh->map(AmanziMesh::FACE, false);
      
      for (int bf=0; bf!=bf_map.NumMyElements(); ++bf) {
        auto f = f_map.LID(bf_map.GID(bf));
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana->pressure_exact(mesh->face_centroid(f), 0.0);
      }
    } else {
      assert(false && "OperatorDiffusion test harness not implemented for this kind.");
    }
  }

  void setBCsDirichletNeumannBox() {
    auto& bc_value = bc->bc_value();
    auto& bc_model = bc->bc_model();
    
    if (bc->kind() == AmanziMesh::FACE) {
      auto bf_map = mesh->map(AmanziMesh::BOUNDARY_FACE, false);
      auto f_map = mesh->map(AmanziMesh::FACE, false);
      
      for (int bf=0; bf!=bf_map.NumMyElements(); ++bf) {
        auto f = f_map.LID(bf_map.GID(bf));
        auto fc = mesh->face_centroid(f);
        if (fc[0] < 1.e-6) {
          bool flag;
          AmanziGeometry::Point normal = FaceNormalExterior(*mesh, f, &flag);
          AmanziGeometry::Point xf = mesh->face_centroid(f);
          bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_value[f] = ana->velocity_exact(xf, 0.0) * normal / mesh->face_area(f);
        } else {
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = ana->pressure_exact(mesh->face_centroid(f), 0.0);
        }          
      }
    } else {
      assert(false && "OperatorDiffusion test harness not implemented for this kind.");
    }
  }

  void setup(const std::string& pc_name_, bool symmetric_) {
    symmetric = symmetric_;
    pc_name = pc_name_;
    global_op = op->global_operator();
    solution = Teuchos::rcp(new CompositeVector(global_op->DomainMap()));
    solution->PutScalar(0.);
    
    // get and assemble the global operator
    if (pc_name != "identity") {
      global_op->SymbolicAssembleMatrix();
    }

    // create preconditoner using the base operator class
    Teuchos::ParameterList slist = plist.sublist("preconditioners").sublist(pc_name);
    global_op->InitializePreconditioner(slist);

    if (symmetric) {
      Teuchos::ParameterList lop_list = plist.sublist("solvers")
                               .sublist("AztecOO CG").sublist("pcg parameters");
      solver = Teuchos::rcp(new AmanziSolvers::LinearOperatorPCG<Operators::Operator, CompositeVector, CompositeVectorSpace>(global_op, global_op));
      solver->Init(lop_list);
    } else {
      Teuchos::ParameterList lop_list = plist.sublist("solvers")
                               .sublist("GMRES").sublist("GMRES parameters");
      solver = Teuchos::rcp(new AmanziSolvers::LinearOperatorGMRES<Operators::Operator, CompositeVector, CompositeVectorSpace>(global_op, global_op));
      solver->Init(lop_list);
    }

    CompositeVectorSpace flux_space;
    flux_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    flux = flux_space.Create();
  }

  void go(double tol=1.e-14) {
    global_op->Init();
    op->UpdateMatrices(Teuchos::null, solution.ptr());

    CompositeVector& rhs = *global_op->rhs();
    Epetra_MultiVector& rhs_c = *rhs.ViewComponent("cell", false);
    for (int c=0; c!=mesh->num_entities(AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_type::OWNED); ++c) {
      const auto& xc = mesh->cell_centroid(c);
      rhs_c[0][c] += ana->source_exact(xc, 0.0) * mesh->cell_volume(c);
    }
    
    op->ApplyBCs(true, true, true);
    if (pc_name != "identity") {
      global_op->SymbolicAssembleMatrix();
      global_op->AssembleMatrix();
    }
    global_op->UpdatePreconditioner();

    if (symmetric) {
      // Test SPD properties of the preconditioner.
      VerificationCV ver(global_op);
      ver.CheckPreconditionerSPD();
    }

    int ierr = solver->ApplyInverse(*global_op->rhs(), *solution);

    auto MyPID = comm->MyPID();
    if (MyPID == 0) {
      std::cout << "pressure solve: ||r||=" << solver->residual() 
                << " itr=" << solver->num_itrs()
                << " code=" << solver->returned_code() << std::endl;
    }
      
    if (tol > 0.0) {
      // compute pressure error
      double pnorm(0.), pl2_err(0.), pinf_err(0.);
      Epetra_MultiVector& p = *solution->ViewComponent("cell", false);
      ComputeCellError(*ana, mesh, p, 0.0, pnorm, pl2_err, pinf_err);

      // calculate flux error
      Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
      double unorm, ul2_err, uinf_err;
      
      op->UpdateFlux(solution.ptr(), flux.ptr());
      ComputeFaceError(*ana, mesh, flx, 0.0, unorm, ul2_err, uinf_err);
      //flux->Print(std::cout);

      if (MyPID == 0) {
        pl2_err /= pnorm; 
        ul2_err /= unorm;
        printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n",
               pl2_err, pinf_err, ul2_err, uinf_err, solver->num_itrs());
        
        CHECK(pl2_err < tol);
        CHECK(ul2_err < 10*tol);
        if (pc_name != "identity" && pc_name != "diagonal") CHECK(solver->num_itrs() < 10);
      }
    }
  }
  
  Comm_ptr_type comm;
  Teuchos::ParameterList plist;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;

  bool symmetric;
  Teuchos::RCP<Operators::BCs> bc;
  Teuchos::RCP<AnalyticBase> ana;
  Teuchos::RCP<CompositeVector> solution, flux;
  std::string pc_name;
  
  Teuchos::RCP<Operators::PDE_Diffusion> op;
  Teuchos::RCP<Operators::Operator> global_op;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace>> solver;
};


#endif
