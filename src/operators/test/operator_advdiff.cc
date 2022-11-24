/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionMFD.hh"

#include "Analytic00.hh"


/* *****************************************************************
* Verify convergence of advection-diffusion solver for various BCs.
* **************************************************************** */
template <class Analytic>
void
AdvectionDiffusionStatic(int nx, double* error)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Advection-duffusion in 2D" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_advdiff.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, nx);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // populate diffusion coefficient
  AmanziGeometry::Point vel(1.0, 1.0);
  Analytic ana(mesh, 2, 0.0, vel);

  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells_owned; c++) {
    const auto& xc = mesh->cell_centroid(c);
    auto Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create discretization space
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);

  // create source
  CompositeVector source(*cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const auto& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0) * mesh->cell_volume(c);
  }

  // create flux field
  auto u = Teuchos::rcp(new CompositeVector(*cvs));
  Epetra_MultiVector& uf = *u->ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    auto tmp = ana.advection_exact(mesh->face_centroid(f), 0.0);
    uf[0][f] = tmp * mesh->face_normal(f);
  }

  // create boundary data for diffusion and
  // Neumann on outflow boundary is needed for operator's positive-definitness
  bool flag;
  double diff_flux, adv_flux;
  auto bcd = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  {
    std::vector<int>& bc_model = bcd->bc_model();
    std::vector<double>& bc_value = bcd->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const auto& xf = mesh->face_centroid(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc_model[f] = OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana.pressure_exact(xf, 0.0);
      } else if (fabs(xf[1]) < 1e-6) { // inflow bottom boundary
        double area = mesh->face_area(f);
        auto normal = ana.face_normal_exterior(f, &flag) / area;

        bc_model[f] = OPERATOR_BC_TOTAL_FLUX;
        diff_flux = ana.velocity_exact(xf, 0.0) * normal;
        adv_flux = (ana.advection_exact(xf, 0.0) * normal) * ana.pressure_exact(xf, 0.0);
        bc_value[f] = diff_flux + adv_flux;
      } else if (fabs(xf[0] - 1.0) < 1e-6) { // outflow right boundary
        double area = mesh->face_area(f);
        auto normal = ana.face_normal_exterior(f, &flag);

        bc_model[f] = OPERATOR_BC_NEUMANN;
        bc_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
      }
    }
  }

  // This split of BCs was designed for testing purposes only!
  // we need to specify only the Dirichlet boundary condition on inflow boundary
  // Dirichlet on the outlow boundary is ignored by upwind operator
  // Neumann condition on inflow violates monotonicity; it generates negative
  //   contribution to diagonal
  auto bca = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  {
    std::vector<int>& bc_model = bca->bc_model();
    std::vector<double>& bc_value = bca->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const auto& xf = mesh->face_centroid(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc_model[f] = OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana.pressure_exact(xf, 0.0);
      } else if (fabs(xf[1]) < 1e-6) { // inflow bottom boundary
        bc_model[f] = OPERATOR_BC_TOTAL_FLUX;
        // bc_value[f] = not used
      } else if (fabs(xf[0] - 1.0) < 1e-6) { // outflow right boundary
        bc_model[f] = OPERATOR_BC_NEUMANN;
        // bc_value[f] = not used
      }
    }
  }
  // loop until convergence
  // create diffusion operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator");
  auto op_diff = Teuchos::rcp(new PDE_DiffusionMFD(olist, mesh));
  op_diff->Init(olist);
  op_diff->SetBCs(bcd, bcd);

  // set up the diffusion operator
  op_diff->Setup(K, Teuchos::null, Teuchos::null);
  op_diff->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_diff->ApplyBCs(true, true, true);

  // create an advection operator
  // this is minor op for advection-mdiffusion pair; hence, BCs are secondary
  auto global_op = op_diff->global_operator();
  olist = plist.sublist("PK operator").sublist("advection operator");
  auto op_adv = Teuchos::rcp(new PDE_AdvectionUpwind(olist, global_op));
  op_adv->SetBCs(bca, bca);

  op_adv->Setup(*u);
  op_adv->UpdateMatrices(u.ptr());
  op_adv->ApplyBCs(false, true, false);

  global_op->set_inverse_parameters(
    "Hypre AMG", plist.sublist("preconditioners"), "NKA", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector& rhs = *global_op->rhs();
  rhs.Update(1.0, source, 1.0);

  global_op->ApplyInverse(rhs, solution);

  // compute pressure error
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
  double pnorm, pl2_err, pinf_err;
  ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);
  *error = pl2_err;

  if (MyPID == 0) {
    pl2_err /= pnorm;
    printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, global_op->num_itrs());

    CHECK(pl2_err < 3e-2);
    CHECK(global_op->num_itrs() < 10);
  }

  // visualization
  if (MyPID == 0) {
    std::cout << "pressure solver (gmres): ||r||=" << global_op->residual()
              << " itr=" << global_op->num_itrs() << " code=" << global_op->returned_code()
              << std::endl;

    // visualization
    const Epetra_MultiVector& pc = *solution.ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string) "operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(pc, 0, "solution");
    GMV::close_data_file();
  }
}

TEST(ADVECTION_DIFFUSION_STATIC)
{
  double err1, err2;
  AdvectionDiffusionStatic<Analytic00>(10, &err1);
  AdvectionDiffusionStatic<Analytic00>(20, &err2);
  CHECK(err2 < err1 / 1.8);
}


/* *****************************************************************
* Verify convergence of advection-diffusion solver for various BCs.
* **************************************************************** */
template <class Analytic>
void
AdvectionDiffusionTransient(int nx, int nloops, double* error)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Advection-duffusion transient" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_advdiff.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, nx);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // populate diffusion coefficient
  AmanziGeometry::Point vel(10.0, 0.0);
  Analytic ana(mesh, 1, 0.0, vel);

  Teuchos::RCP<std::vector<WhetStone::Tensor>> K =
    Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  for (int c = 0; c < ncells_owned; c++) {
    const auto& xc = mesh->cell_centroid(c);
    auto Kc = ana.TensorDiffusivity(xc, 0.0);
    K->push_back(Kc);
  }

  // create discretization space
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  CompositeVector solution(*cvs);
  solution.PutScalar(0.0);

  // create source
  CompositeVector source(*cvs);
  Epetra_MultiVector& src = *source.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const auto& xc = mesh->cell_centroid(c);
    src[0][c] = ana.source_exact(xc, 0.0) * mesh->cell_volume(c);
  }

  // create flux field
  auto cvs_u = Teuchos::rcp(new CompositeVectorSpace());
  cvs_u->SetMesh(mesh)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);
  auto u = Teuchos::rcp(new CompositeVector(*cvs_u));
  Epetra_MultiVector& uf = *u->ViewComponent("face");
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    auto tmp = ana.advection_exact(mesh->face_centroid(f), 0.0);
    uf[0][f] = tmp * mesh->face_normal(f);
  }

  // create boundary data for diffusion and
  // Neumann on outflow boundary is needed for operator's positive-definitness
  bool flag;
  auto bcd = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  {
    std::vector<int>& bc_model = bcd->bc_model();
    std::vector<double>& bc_value = bcd->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const auto& xf = mesh->face_centroid(f);
      if (fabs(xf[0]) < 1e-6) {
        bc_model[f] = OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana.pressure_exact(xf, 0.0);
      }
      // else {
      else if (fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        double area = mesh->face_area(f);
        auto normal = ana.face_normal_exterior(f, &flag);

        bc_model[f] = OPERATOR_BC_NEUMANN;
        bc_value[f] = (ana.velocity_exact(xf, 0.0) * normal) / area;
      }
    }
  }

  // This split of BCs was designed for testing purposes only!
  // we need to specify only the Dirichlet boundary condition on inflow boundary
  // Dirichlet on the outlow boundary is ignored by upwind operator
  // Neumann condition on inflow violates monotonicity; it generates negative
  //   contribution to diagonal
  auto bca = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  {
    std::vector<int>& bc_model = bca->bc_model();
    std::vector<double>& bc_value = bca->bc_value();

    for (int f = 0; f < nfaces_wghost; f++) {
      const auto& xf = mesh->face_centroid(f);
      if (fabs(xf[0]) < 1e-6) {
        bc_model[f] = OPERATOR_BC_DIRICHLET;
        bc_value[f] = ana.pressure_exact(xf, 0.0);
      }
      // else {
      else if (fabs(xf[0] - 1.0) < 1e-6 || // outflow right boundary
               fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc_model[f] = OPERATOR_BC_NEUMANN;
        // bc_value[f] = not used
      }
    }
  }

  // loop until convergence
  Teuchos::RCP<Operator> global_op;
  for (int loop = 0; loop < nloops; ++loop) {
    // create diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("diffusion operator fv");
    auto op_diff = Teuchos::rcp(new PDE_DiffusionFV(olist, mesh));
    op_diff->SetBCs(bcd, bcd);

    // set up the diffusion operator
    op_diff->Setup(K, Teuchos::null, Teuchos::null);
    op_diff->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_diff->ApplyBCs(true, true, true);

    // create an advection operator
    // this is minor op for advection-mdiffusion pair; hence, BCs are secondary
    global_op = op_diff->global_operator();
    olist = plist.sublist("PK operator").sublist("advection operator");
    auto op_adv = Teuchos::rcp(new PDE_AdvectionUpwind(olist, global_op));
    op_adv->SetBCs(bca, bca);

    op_adv->Setup(*u);
    op_adv->UpdateMatrices(u.ptr());
    op_adv->ApplyBCs(false, true, false);

    // add accumulation terms
    if (nloops > 1) {
      double dT(1.0);
      CompositeVector phi(*cvs);
      phi.PutScalar(0.2);

      PDE_Accumulation op_acc(AmanziMesh::CELL, global_op);
      op_acc.AddAccumulationDelta(solution, phi, phi, dT, "cell");
    }

    global_op->set_inverse_parameters(
      "Hypre AMG", plist.sublist("preconditioners"), "AztecOO CG", plist.sublist("solvers"));
    global_op->InitializeInverse();
    global_op->ComputeInverse();

    CompositeVector& rhs = *global_op->rhs();
    rhs.Update(1.0, source, 1.0);

    global_op->ApplyInverse(rhs, solution);

    // compute pressure error
    Epetra_MultiVector& p = *solution.ViewComponent("cell", false);
    double pnorm, pl2_err, pinf_err;
    ana.ComputeCellError(p, 0.0, pnorm, pl2_err, pinf_err);
    *error = pl2_err;

    if (MyPID == 0) {
      pl2_err /= pnorm;
      printf("L2(p)=%9.6f  Inf(p)=%9.6f  itr=%3d\n", pl2_err, pinf_err, global_op->num_itrs());

      if (loop == nloops - 1) {
        CHECK(pl2_err < 3e-2);
        CHECK(global_op->num_itrs() < 10);
      }
    }
  }
}


TEST(ADVECTION_DIFFUSION_TRANSIENT)
{
  double err1;
  AdvectionDiffusionTransient<Analytic00>(40, 5, &err1);
}
