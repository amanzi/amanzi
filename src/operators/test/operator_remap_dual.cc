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
#include "MeshMaps_FEM.hh"
#include "MeshMaps_VEM.hh"
#include "Tensor.hh"
#include "WhetStone_typedefs.hh"

// Amanzi::Operators
#include "Abstract.hh"
#include "Accumulation.hh"
#include "AdvectionRiemann.hh"
#include "OperatorDefs.hh"
#include "Reaction.hh"
#include "RemapUtils.hh"


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapTests2DDual(int dim, int order, std::string disc_name,
                      std::string maps_name,
                      int nx, int ny, double dt) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: " << dim << "D remap, dual formulation: order=" 
                            << order << ", maps=" << maps_name << std::endl;

  // polynomial space
  WhetStone::Polynomial pp(dim, order);
  int nk = pp.size();

  // create initial mesh
  MeshFactory meshfactory(&comm);
  meshfactory.set_partitioner(AmanziMesh::Partitioner_type::METIS);
  meshfactory.preference(FrameworkPreference({MSTK}));

  Teuchos::RCP<const Mesh> mesh0;
  if (dim == 2) {
    // mesh0 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
    mesh0 = meshfactory("test/median15x16.exo", Teuchos::null);
  } else {
    mesh0 = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, ny, Teuchos::null, true, true);
  }

  int ncells_owned = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh0->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int nfaces_owned = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh0->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nnodes_owned = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh0->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  // create second and auxiliary mesh
  Teuchos::RCP<Mesh> mesh1;
  if (dim == 2) {
    // mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
    mesh1 = meshfactory("test/median15x16.exo", Teuchos::null);
  } else {
    mesh1 = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, ny, Teuchos::null, true, true);
  }

  // deform the second mesh
  AmanziGeometry::Point xv(dim), xref(dim), uv(dim);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh1->node_get_coordinates(v, &xv);

    double ds(0.001);
    for (int i = 0; i < 1000; ++i) {
      if (dim == 2) {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
      } else {
        // uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
        // uv[1] =-0.1 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
        // uv[2] =-0.1 * std::cos(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::sin(M_PI * xv[2]);
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
        uv[2] = 0.0;
      }

      xv += uv * ds;
    }

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh1->deform(nodeids, new_positions, false, &final_positions);

  // little factory of mesh maps
  std::shared_ptr<WhetStone::MeshMaps> maps;
  if (maps_name == "FEM") {
    maps = std::make_shared<WhetStone::MeshMaps_FEM>(mesh0, mesh1);
  } else if (maps_name == "VEM") {
    maps = std::make_shared<WhetStone::MeshMaps_VEM>(mesh0, mesh1);
  }

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  // we need dg to compute scaling of basis functions
  WhetStone::DG_Modal dg(order, mesh0);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh0->cell_centroid(c);
    p1c[0][c] = std::sin(3 * xc[0]) * std::sin(6 * xc[1]);
    if (nk > 1) {
      double a, b;
      WhetStone::Iterator it(dim);

      it.begin(1);
      dg.TaylorBasis(c, it, &a, &b);
      p1c[1][c] = 3 * std::cos(3 * xc[0]) * std::sin(6 * xc[1]) / a;

      ++it;
      dg.TaylorBasis(c, it, &a, &b);
      p1c[2][c] = 6 * std::sin(3 * xc[0]) * std::cos(6 * xc[1]) / a;
    }
if (dim == 3) p1c[0][c] = 1.0;
  }

  // initial mass
  double mass0(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    mass0 += p1c[0][c] * mesh0->cell_volume(c);
  }
  double mass_tmp(mass0);
  mesh0->get_comm()->SumAll(&mass_tmp, &mass0, 1);

  // allocate memory
  CompositeVectorSpace cvs2;
  cvs2.SetMesh(mesh1)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  CompositeVector p2(cvs2);
  Epetra_MultiVector& p2c = *p2.ViewComponent("cell");

  // create flux operator
  Teuchos::ParameterList plist;
  plist.set<std::string>("method", disc_name)
       .set<int>("method order", order)
       .set<bool>("jump operator on test function", true);

  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op = Teuchos::rcp(new AdvectionRiemann(plist, mesh0));
  auto global_op = op->global_operator();

  std::vector<WhetStone::VectorPolynomial> vec_vel(nfaces_wghost);
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > vel = 
      Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));

  // Attach volumetric advection operator to the flux operator.
  // We modify the existing parameter list.
  plist.set<std::string>("matrix type", "advection")
       .set<bool>("gradient operator on test function", true);
  plist.sublist("schema domain").set<std::string>("base", "cell");
  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<Abstract> op_adv = Teuchos::rcp(new Abstract(plist, global_op));

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > cell_vel = 
      Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned));
 
  op_adv->SetupPolyVector(cell_vel);

  // approximate accumulation term using the reaction operator
  // -- left-hand-side
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<Reaction> op_reac0 = Teuchos::rcp(new Reaction(plist, mesh0));
  auto global_reac0 = op_reac0->global_operator();

  Teuchos::RCP<WhetStone::VectorPolynomial> jac0 = 
     Teuchos::rcp(new WhetStone::VectorPolynomial(ncells_owned));

  op_reac0->Setup(jac0);

  // -- right-hand-side
  Teuchos::RCP<Reaction> op_reac1 = Teuchos::rcp(new Reaction(plist, mesh0));
  auto global_reac1 = op_reac1->global_operator();

  Teuchos::RCP<WhetStone::VectorPolynomial> jac1 = 
     Teuchos::rcp(new WhetStone::VectorPolynomial(ncells_owned));

  op_reac1->Setup(jac1);

  // explicit time integration
  double t(0.0), tend(1.0);
  while(t < tend - dt/2) {
    // calculate face velocities
    for (int f = 0; f < nfaces_wghost; ++f) {
      maps->VelocityFace(f, vec_vel[f]);
    }

    // calculate normal component of face velocity at time t+dt/2 
    for (int f = 0; f < nfaces_wghost; ++f) {
      // cn = j J^{-t} N dA
      WhetStone::VectorPolynomial cn;
      maps->NansonFormula(f, t + dt/2, vec_vel[f], cn);
      (*vel)[f] = vec_vel[f] * cn;
    }

    // calculate cell velocity at time t+dt/2 and change its sign
    Entity_ID_List faces;
    WhetStone::MatrixPolynomial C;
    WhetStone::VectorPolynomial tmp;

    for (int c = 0; c < ncells_owned; ++c) {
      mesh0->cell_get_faces(c, &faces);
      std::vector<WhetStone::VectorPolynomial> vf;

      for (int n = 0; n < faces.size(); ++n) {
        vf.push_back(vec_vel[faces[n]]);
      }

      maps->VelocityCell(c, vf, tmp);
      maps->Cofactors(c, t + dt/2, tmp, C);
      tmp[0].Multiply(C, tmp, (*cell_vel)[c], true);

      for (int i = 0; i < dim; ++i) {
        (*cell_vel)[c][i] *= -1.0;
      }
    }

    // calculate determinant of Jacobian at time t+dt
    for (int c = 0; c < ncells_owned; ++c) {
      mesh0->cell_get_faces(c, &faces);
      std::vector<WhetStone::VectorPolynomial> vf;

      for (int n = 0; n < faces.size(); ++n) {
        vf.push_back(vec_vel[faces[n]]);
      }

      maps->JacobianDet(c, t, vf, (*jac0)[c]);
      maps->JacobianDet(c, t + dt, vf, (*jac1)[c]);
    }

    // populate operators
    op->UpdateMatrices(vel.ptr());
    op_adv->UpdateMatrices();
    op_reac0->UpdateMatrices(p1.ptr());
    op_reac1->UpdateMatrices(p1.ptr());

    // predictor step
    // -- calculate rhs
    CompositeVector& rhs = *global_reac0->rhs();
    global_reac0->Apply(*p1, rhs);

    CompositeVector g(cvs1);
    global_op->Apply(*p1, g);
    g.Update(1.0, rhs, dt);

    // -- solve the problem with mass matrix
    global_reac1->SymbolicAssembleMatrix();
    global_reac1->AssembleMatrix();

    plist.set<std::string>("preconditioner type", "diagonal");
    global_reac1->InitPreconditioner(plist);

    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        pcg(global_reac1, global_reac1);

    std::vector<std::string> criteria;
    criteria.push_back("absolute residual");
    plist.set<double>("error tolerance", 1e-12)
         .set<Teuchos::Array<std::string> >("convergence criteria", criteria);
    pcg.Init(plist);
    pcg.ApplyInverse(g, p2);

    // corrector step
    p2.Update(0.5, *p1, 0.5);
    global_op->Apply(p2, g);
    g.Update(1.0, rhs, dt);

    pcg.ApplyInverse(g, p2);

    // end timestep operations
    *p1->ViewComponent("cell") = *p2.ViewComponent("cell");
    t += dt;
  }

  // calculate error in the new basis
  double mass1(0.0);
  double pl2_err(0.0), pinf_err(0.0), area(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    double area_c(mesh1->cell_volume(c));
    if (nk == 1) {
      // const AmanziGeometry::Point& xg = maps->cell_geometric_center(1, c);
      const AmanziGeometry::Point& xg = mesh1->cell_centroid(c);
      double tmp = p2c[0][c] - std::sin(3 * xg[0]) * std::sin(6 * xg[1]);
if (dim == 3) tmp = p2c[0][c] - 1.0;

      pinf_err = std::max(pinf_err, fabs(tmp));
      pl2_err += tmp * tmp * area_c;
    }
    else {
      std::vector<double> data;
      for (int i = 0; i < nk; ++i) {
        data.push_back(p2c[i][c]);
      }
      WhetStone::Polynomial poly(dg.CalculatePolynomial(c, data));

      Entity_ID_List nodes;
      AmanziGeometry::Point v0(dim), v1(dim);

      mesh0->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();  
      for (int i = 0; i < nnodes; ++i) {
        mesh0->node_get_coordinates(nodes[i], &v0);
        mesh1->node_get_coordinates(nodes[i], &v1);

        double tmp = poly.Value(v0);
        tmp -= std::sin(3 * v1[0]) * std::sin(6 * v1[1]);
        pinf_err = std::max(pinf_err, fabs(tmp));
        pl2_err += tmp * tmp * area_c / nnodes;
      }
    }

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
  CHECK(pl2_err < 0.12 / (order + 1));

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

/*
TEST(REMAP2D_DG0_DUAL_FEM) {
  RemapTests2DDual(2, 0, "dg modal", "FEM", 10, 10, 0.1);
}

TEST(REMAP2D_DG1_DUAL_FEM) {
  RemapTests2DDual(2, 1, "dg modal", "FEM", 10, 10, 0.1);
}
*/

TEST(REMAP2D_DG0_DUAL_VEM) {
  RemapTests2DDual(2, 0, "dg modal", "VEM", 10, 10, 0.05);
}

TEST(REMAP2D_DG1_DUAL_VEM) {
  RemapTests2DDual(2, 1, "dg modal", "VEM", 10, 10, 0.05);
}
/*

TEST(REMAP3D_DG0_DUAL_VEM) {
  RemapTests2DDual(3, 0, "dg modal", "VEM", 10, 10, 0.1);
}
*/

