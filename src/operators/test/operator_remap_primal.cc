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
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "DG_Modal.hh"
#include "GMVMesh.hh"
#include "LinearOperatorPCG.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "MeshUtils.hh"
#include "Tensor.hh"
#include "WhetStone_typedefs.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

#include "AnalyticDG04.hh"

/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Primal formulation uses gradient and jumps of a solution.
***************************************************************** */
void RemapTests2DPrimal(int order, std::string disc_name,
                        std::string maps_name,
                        int nx, int ny, double dt) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: remap in 2D, order=" << order 
                            << ", maps=" << maps_name << std::endl;

  // polynomial space
  int nk = (order + 1) * (order + 2) / 2;

  // create initial mesh
  MeshFactory meshfactory(&comm);
  meshfactory.set_partitioner(AmanziMesh::Partitioner_type::METIS);
  meshfactory.preference(FrameworkPreference({MSTK}));

  Teuchos::RCP<const Mesh> mesh0 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
  // Teuchos::RCP<const Mesh> mesh0 = meshfactory("test/median15x16.exo", Teuchos::null);

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_owned = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nnodes_wghost = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // create second and auxiliary mesh
  Teuchos::RCP<Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
  // Teuchos::RCP<Mesh> mesh1 = meshfactory("test/median15x16.exo", Teuchos::null);

  // deform the second mesh
  AmanziGeometry::Point xv(2), xref(2);
  Entity_ID_List nodeids, faces;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    double ds(0.001), ux, uy;
    for (int i = 0; i < 1000; ++i) {
      ux = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
      uy =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);

      xv[0] += ux * ds;
      xv[1] += uy * ds;
    }

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // little factory of mesh maps
  Teuchos::ParameterList map_list;
  map_list.set<std::string>("method", "CrouzeixRaviart")
          .set<int>("method order", order + 1)
          .set<std::string>("projector", "H1 harmonic")
          .set<std::string>("map name", maps_name);

  WhetStone::MeshMapsFactory maps_factory;
  auto maps = maps_factory.Create(map_list, mesh0, mesh1);

  // numerical integration
  WhetStone::NumericalIntegration numi(mesh0);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  // we need dg to compute scaling of basis functions
  WhetStone::DG_Modal dg(order, mesh0, "orthonormalized");

  AnalyticDG04 ana(mesh0, order);
  ana.InitialGuess(dg, p1c, 1.0);

  // initial mass
  double mass0(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    mass0 += p1c[0][c] * mesh0->cell_volume(c);
  }
  double mass_tmp(mass0);
  mesh0->get_comm()->SumAll(&mass_tmp, &mass0, 1);

  // allocate memory for global variables
  // -- solution
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // -- face velocities
  std::vector<WhetStone::VectorPolynomial> vec_vel(nfaces_wghost);
  for (int f = 0; f < nfaces_wghost; ++f) {
    maps->VelocityFace(f, vec_vel[f]);
  }

  // -- cell-baced velocities and Jacobian matrices
  std::vector<WhetStone::VectorPolynomial> uc(ncells_owned);
  std::vector<WhetStone::MatrixPolynomial> J(ncells_owned);

  for (int c = 0; c < ncells_owned; ++c) {
    mesh0->cell_get_faces(c, &faces);

    std::vector<WhetStone::VectorPolynomial> vvf;
    for (int n = 0; n < faces.size(); ++n) {
      vvf.push_back(vec_vel[faces[n]]);
    }

    maps->VelocityCell(c, vvf, uc[c]);
    maps->Jacobian(uc[c], J[c]);
  }

  // create flux operator
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", disc_name)
       .set<std::string>("dg basis", "orthonormalized")
       .set<std::string>("matrix type", "flux")
       .set<std::string>("flux formula", "downwind")
       .set<int>("method order", order)
       .set<bool>("jump operator on test function", false);

  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<PDE_AdvectionRiemann> op = Teuchos::rcp(new PDE_AdvectionRiemann(plist, mesh0));
  auto global_op = op->global_operator();

  Teuchos::RCP<std::vector<WhetStone::Polynomial> > vel = 
      Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));

  // Attach volumetric advection operator to the flux operator.
  // We modify the existing parameter list.
  plist.set<std::string>("matrix type", "advection")
       .set<bool>("gradient operator on test function", false);
  plist.sublist("schema domain").set<std::string>("base", "cell");
  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<PDE_Abstract> op_adv = Teuchos::rcp(new PDE_Abstract(plist, global_op));

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > cell_vel = 
      Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned));
 
  op_adv->SetupPolyVector(cell_vel);

  // create accumulation operator re-using reaction operator
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<PDE_Reaction> op_reac = Teuchos::rcp(new PDE_Reaction(plist, mesh0));
  auto global_reac = op_reac->global_operator();

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > jac = 
     Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned));

  op_reac->Setup(jac);

  double t(0.0), tend(1.0);
  while(t < tend - dt/2) {
    // calculate normal component of face velocity at time t+dt/2 
    for (int f = 0; f < nfaces_wghost; ++f) {
      // cn = j J^{-t} N dA
      WhetStone::VectorPolynomial cn;
      maps->NansonFormula(f, t + dt/2, vec_vel[f], cn);

      (*vel)[f] = vec_vel[f] * cn;
      (*vel)[f] *= -1.0;
    }

    // calculate various geometric quantities in reference framework
    WhetStone::MatrixPolynomial C;

    for (int c = 0; c < ncells_owned; ++c) {
      maps->Cofactors(t + dt/2, J[c], C);
      (*cell_vel)[c].Multiply(C, uc[c], true);

      maps->Determinant(t + dt, J[c], (*jac)[c]);
    }

    // populate operators
    op->UpdateMatrices(vel.ptr());
    op_adv->UpdateMatrices();
    op_reac->UpdateMatrices(p1.ptr());

    // predictor step
    // -- calculate rhs
    CompositeVector& rhs = *global_reac->rhs();
    global_reac->Apply(*p1, rhs);

    CompositeVector g(cvs1);
    global_op->Apply(*p1, g);
    g.Update(1.0, rhs, dt);

    // -- solve the problem with mass matrix
    global_reac->SymbolicAssembleMatrix();
    global_reac->AssembleMatrix();

    plist.set<std::string>("preconditioner type", "diagonal");
    global_reac->InitializePreconditioner(plist);
    global_reac->UpdatePreconditioner();

    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        pcg(global_reac, global_reac);

    plist.set<double>("error tolerance", 1e-12);
    pcg.Init(plist);
    pcg.ApplyInverse(g, p2);

    // corrector step
    /*
    p2.Update(0.5, p1, 0.5);
    global_op->Apply(p2, g);
    g.Update(1.0, rhs, 1.0);

    pcg.ApplyInverse(g, p2);
    */

    *p1->ViewComponent("cell") = *p2.ViewComponent("cell");
    t += dt;
  }

  // calculate error in the new basis
  double mass1(0.0);
  double pl2_err(0.0), pinf_err(0.0), area(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    AmanziGeometry::Point xg = WhetStone::cell_geometric_center(*mesh1, c);
    double area_c = mesh1->cell_volume(c);

    double tmp = p2c[0][c] - std::sin(3 * xg[0]) * std::sin(6 * xg[1]);
    pinf_err = std::max(pinf_err, fabs(tmp));
    pl2_err += tmp * tmp * area_c;

    area += area_c;
    mass1 += p2c[0][c] * mesh1->cell_volume(c);
  }

  // parallel colelctive operations
  double err_tmp(pl2_err);
  mesh1->get_comm()->SumAll(&err_tmp, &pl2_err, 1);

  err_tmp = area;
  mesh1->get_comm()->SumAll(&err_tmp, &area, 1);

  err_tmp = pinf_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &pinf_err, 1);

  mass_tmp = mass1;
  mesh1->get_comm()->SumAll(&mass_tmp, &mass1, 1);

  // error tests
  pl2_err = std::pow(pl2_err, 0.5);
  if (MyPID == 0) CHECK(pl2_err < 0.08 / (order + 1));

  if (MyPID == 0) {
    printf("L2(p0)=%12.8g  Inf(p0)=%12.8g  dMass=%12.8g  Err(area)=%12.8g\n", 
        pl2_err, pinf_err, mass1 - mass0, 1.0 - area);
  }

  // visualization
  if (MyPID == 0) {
    GMV::open_data_file(*mesh1, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p2c, 0, "remaped");
    if (order > 0) {
      GMV::write_cell_data(p2c, 1, "gradx");
      GMV::write_cell_data(p2c, 2, "grady");
    }
    GMV::close_data_file();
  }
}


TEST(REMAP_DG0_PRIMAL_FEM) {
  RemapTests2DPrimal(0, "dg modal", "FEM", 12, 12, 0.1);
}

TEST(REMAP_DG1_PRIMAL_FEM) {
  RemapTests2DPrimal(1, "dg modal", "FEM", 10, 10, 0.1);
}

TEST(REMAP_DG0_PRIMAL_VEM) {
  RemapTests2DPrimal(0, "dg modal", "VEM", 10, 10, 0.1);
}

TEST(REMAP_DG1_PRIMAL_VEM) {
  RemapTests2DPrimal(1, "dg modal", "VEM", 10, 10, 0.1);
}


