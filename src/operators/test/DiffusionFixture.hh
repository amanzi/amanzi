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
#include "IterativeMethodGMRES.hh"
#include "IterativeMethodNKA.hh"
#include "IterativeMethodPCG.hh"

#include "Tensor.hh"


#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"

using namespace Amanzi;

struct DiffusionFixture {
  DiffusionFixture(const Teuchos::RCP<AnalyticBase>& ana_,
                   const std::string& mesh_file,
                   Teuchos::RCP<const AmanziMesh::Mesh> mesh_ = Teuchos::null)
      : comm(Amanzi::getDefaultComm()),
        ana(ana_),
        mesh(mesh_)
  {
    plist = Teuchos::getParametersFromXmlFile("test/operator_diffusion_discretizations.xml");

    if (mesh == Teuchos::null ){
      auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(ana->dimension(),
            plist->sublist("regions"), *comm));
      AmanziMesh::MeshFactory meshfactory(comm, gm);
      if (mesh_file == "Generate1D") {
        mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 100, 1);
      } else if (mesh_file == "Generate2D") {
        mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 100, 10);
      } else if (mesh_file == "Generate2D_HiRes") {
        mesh = meshfactory.create(-1.0, -1.0, 1.0, 1.0, 100, 100);
      } else if (mesh_file == "Generate3D") {
        mesh = meshfactory.create(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 4, 5, 6);
      } else {
        mesh = meshfactory.create(mesh_file);
      }
    }

  }

  template<class PDE_Diffusion_type, AmanziMesh::Entity_kind Boundary_kind>
  void discretize(const std::string& name) {
    // construct the operator
    op = Teuchos::rcp(new PDE_Diffusion_type(plist->sublist("PK operator").sublist(name), mesh));
    op->Init();

    // modify diffusion coefficient
    CompositeVectorSpace K_map;
    K_map.SetMesh(mesh);
    K_map.AddComponent("cell", AmanziMesh::CELL, 1);
    auto K = Teuchos::rcp(new TensorVector(K_map));



    std::vector<WhetStone::Tensor<DefaultHostMemorySpace>> host_tensors(K->size());
    for (int c = 0; c < K->size(); c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      host_tensors[c].assign(ana->TensorDiffusivity(xc, 0.0));
    }
    K->Init(host_tensors);
    op->SetTensorCoefficient(K);
    // boundary condition
    bc = Teuchos::rcp(new Operators::BCs(mesh, Boundary_kind, WhetStone::DOF_Type::SCALAR));
    op->SetBCs(bc,bc);
  }

  template<class PDE_Diffusion_type, AmanziMesh::Entity_kind Boundary_kind>
  void discretizeWithGravity(const std::string& name, double gravity) {
    AmanziGeometry::Point g(mesh->space_dimension());
    g[mesh->space_dimension()-1] = -gravity;
    Teuchos::ParameterList op_list = plist->sublist("PK operator").sublist(name);
    op_list.set("gravity", true);
    op = Teuchos::rcp(new PDE_Diffusion_type(op_list, mesh));
    op->SetDensity(1.0);
    op->SetGravity(g);
    op->Init();

    // modify diffusion coefficient
    CompositeVectorSpace K_map;
    K_map.SetMesh(mesh);
    K_map.AddComponent("cell", AmanziMesh::CELL, 1);
    auto K = Teuchos::rcp(new TensorVector(K_map));

    std::vector<WhetStone::Tensor<DefaultHostMemorySpace>> host_tensors(K->size());
    for (int c = 0; c < K->size(); c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      host_tensors[c].assign(ana->TensorDiffusivity(xc, 0.0));
    }
    K->Init(host_tensors);

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
      auto vec = kr->ViewComponent<MirrorHost>("cell", true);
      for (int i=0; i!=mesh->num_entities(kind, AmanziMesh::Parallel_type::ALL); ++i) {
        vec(i,0) = ana->ScalarDiffusivity(mesh->cell_centroid(i), 0.0);
      }

    } else if (kind == AmanziMesh::FACE) {
      space.SetComponent("face", AmanziMesh::FACE, 1);
      kr = space.Create();
      auto vec = kr->ViewComponent<MirrorHost>("face", true);
      for (int i=0; i!=mesh->num_entities(kind, AmanziMesh::Parallel_type::ALL); ++i) {
        vec(i,0) = ana->ScalarDiffusivity(mesh->face_centroid(i), 0.0);
      }
    }
    op->SetScalarCoefficient(kr, Teuchos::null);
  }

  void setBCsDirichlet() {
    auto bc_value = bc->bc_value();
    auto bc_model = bc->bc_model();

    if (bc->kind() == AmanziMesh::FACE) {
      const auto& bf_map = *mesh->map(AmanziMesh::BOUNDARY_FACE, false);
      const auto& f_map = *mesh->map(AmanziMesh::FACE, false);

      for (int bf=0; bf!=bf_map.getNodeNumElements(); ++bf) {
        auto f = f_map.getLocalElement(bf_map.getGlobalElement(bf));
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana->pressure_exact(mesh->face_centroid(f), 0.0);
      }
    } else {
      assert(false && "OperatorDiffusion test harness not implemented for this kind.");
    }
  }

  void setBCsDirichletNeumannBox() {
    auto bc_value = bc->bc_value();
    auto bc_model = bc->bc_model();

    if (bc->kind() == AmanziMesh::FACE) {
      const auto& bf_map = *mesh->map(AmanziMesh::BOUNDARY_FACE, false);
      const auto& f_map = *mesh->map(AmanziMesh::FACE, false);

      for (int bf=0; bf!=bf_map.getNodeNumElements(); ++bf) {
        auto f = f_map.getLocalElement(bf_map.getGlobalElement(bf));
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

  void setup(const std::string& pc_name_, const std::string& solver_name_) {
    pc_name = pc_name_;
    solver_name = solver_name_;
    global_op = op->global_operator();
    solution = Teuchos::rcp(new CompositeVector(global_op->getDomainMap()));
    solution->putScalar(0.);

    // create preconditoner using the base operator class
    if (pc_name != "none") {
      inverse_list.setParameters(*Teuchos::sublist(Teuchos::sublist(plist, "preconditioners"), pc_name));
    }
    if (solver_name != "none") {
      inverse_list.setParameters(*Teuchos::sublist(Teuchos::sublist(plist, "solvers"), solver_name));
    }
    global_op->set_inverse_parameters(inverse_list);
    global_op->initializeInverse();

    CompositeVectorSpace flux_space;
    flux_space.SetMesh(mesh)->SetGhosted(true)->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    flux = flux_space.Create();
  }

  void go(double tol=1.e-14) {
    global_op->Zero();
    op->UpdateMatrices(Teuchos::null, solution.ptr());

    CompositeVector& rhs = *global_op->rhs();
    auto rhs_c = rhs.ViewComponent<MirrorHost>("cell", false);
    for (int c=0; c!=mesh->num_entities(AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_type::OWNED); ++c) {
      const auto& xc = mesh->cell_centroid(c);
      rhs_c(c,0) += ana->source_exact(xc, 0.0) * mesh->cell_volume(c);
    }
    op->ApplyBCs(true, true, true);

    global_op->computeInverse();
    int ierr = global_op->applyInverse(*global_op->rhs(), *solution);
    //    inverse_list.print(std::cout);

    auto MyPID = comm->getRank();
    if (MyPID == 0) {
      std::cout << "pressure solve: ||r||=" << global_op->residual()
                << " itr=" << global_op->num_itrs()
                << " code=" << global_op->returned_code() << std::endl;
    }

    if (tol > 0.0) {
      // compute pressure error
      double pnorm(0.), pl2_err(0.), pinf_err(0.);
      ComputeCellError(*ana, mesh, *solution, 0.0, pnorm, pl2_err, pinf_err);

      // calculate flux error
      double unorm, ul2_err, uinf_err;

      op->UpdateMatrices(Teuchos::null, Teuchos::null);
      op->UpdateFlux(solution.ptr(), flux.ptr());
      {
        auto fv = flux->ViewComponent<MirrorHost>("face", false);
      }


      ComputeFaceError(*ana, mesh, *flux, 0.0, unorm, ul2_err, uinf_err);
      //flux->Print(std::cout);

      if (MyPID == 0) {
        pl2_err /= pnorm;
        ul2_err /= unorm;
        printf("L2(p)=%9.6f  Inf(p)=%9.6f  L2(u)=%9.6g  Inf(u)=%9.6f itr=%3d\n",
               pl2_err, pinf_err, ul2_err, uinf_err, global_op->num_itrs());

        CHECK(pl2_err < tol);
        CHECK(ul2_err < 10*tol);
        //if (pc_name != "identity" && pc_name != "diagonal") CHECK(global_op->num_itrs() < 10);
      }
    }
  }

  Comm_ptr_type comm;
  ParameterList_ptr_type plist;
  Teuchos::ParameterList inverse_list;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;

  Teuchos::RCP<Operators::BCs> bc;
  Teuchos::RCP<AnalyticBase> ana;
  Teuchos::RCP<CompositeVector> solution, flux;
  std::string pc_name, solver_name;

  Teuchos::RCP<Operators::PDE_Diffusion> op;
  Teuchos::RCP<Operators::Operator> global_op;
  Teuchos::RCP<AmanziSolvers::InverseIterativeMethod<Operators::Operator,Operators::Operator, CompositeVector, CompositeSpace>> solver;
};


#endif
