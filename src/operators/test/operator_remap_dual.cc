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
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"
#include "RemapUtils.hh"

// #define DEBUG

/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapTestsDual(int dim, int order, std::string disc_name,
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
    if (nx != 0) 
      mesh0 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
    else 
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
    if (nx != 0) 
      mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, nx, ny);
    else 
      mesh1 = meshfactory("test/median15x16.exo", Teuchos::null);
  } else {
    mesh1 = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, ny, Teuchos::null, true, true);
  }

  // deform the second mesh
  AmanziGeometry::Point xv(dim), yv(dim), xref(dim), uv(dim);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh1->node_get_coordinates(v, &xv);
    yv = xv;

    double ds(0.001);
    for (int i = 0; i < 1000; ++i) {
      if (dim == 2) {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
      } else {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
        uv[1] =-0.1 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]) * std::cos(M_PI * xv[2]);
        uv[2] =-0.1 * std::cos(M_PI * xv[0]) * std::cos(M_PI * xv[1]) * std::sin(M_PI * xv[2]);
        // uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        // uv[1] =-0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
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

      if (dim == 3) {
        ++it;
        dg.TaylorBasis(c, it, &a, &b);
      }
    }
  }

#ifdef DEBUG
  Teuchos::RCP<CompositeVector> p1_dbg = Teuchos::rcp(new CompositeVector(cvs1));
  *p1_dbg = *p1;
#endif

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

  Teuchos::RCP<PDE_AdvectionRiemann> op = Teuchos::rcp(new PDE_AdvectionRiemann(plist, mesh0));
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

  Teuchos::RCP<PDE_Abstract> op_adv = Teuchos::rcp(new PDE_Abstract(plist, global_op));

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > cell_vel = 
      Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned));
 
  op_adv->SetupPolyVector(cell_vel);

#ifdef DEBUG
  Teuchos::RCP<PDE_Abstract> op_adv_dbg0 = Teuchos::rcp(new PDE_Abstract(plist, mesh0));
  op_adv_dbg0->SetupPolyVector(cell_vel);

  Teuchos::RCP<PDE_Abstract> op_adv_dbg1 = Teuchos::rcp(new PDE_Abstract(plist, mesh0));
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > cell_vel_dbg = 
      Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned));
  op_adv_dbg1->SetupPolyVector(cell_vel_dbg);

  double ame_err(0.0);
#endif

  // approximate accumulation term using the reaction operator
  // -- left-hand-side
  plist.sublist("schema")
      .set<std::string>("base", "cell")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({nk}));

  Teuchos::RCP<PDE_Reaction> op_reac0 = Teuchos::rcp(new PDE_Reaction(plist, mesh0));
  auto global_reac0 = op_reac0->global_operator();

  Teuchos::RCP<std::vector<WhetStone::Polynomial> > jac0 = 
     Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_owned));

  op_reac0->Setup(jac0);

  // -- right-hand-side
  Teuchos::RCP<PDE_Reaction> op_reac1 = Teuchos::rcp(new PDE_Reaction(plist, mesh0));
  auto global_reac1 = op_reac1->global_operator();

  Teuchos::RCP<std::vector<WhetStone::Polynomial> > jac1 = 
     Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_owned));

  op_reac1->Setup(jac1);

#ifdef DEBUG
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > jac0_dbg = 
     Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_owned));
  Teuchos::RCP<PDE_Reaction> op_reac0_dbg = Teuchos::rcp(new PDE_Reaction(plist, mesh0));
  op_reac0_dbg->Setup(jac0_dbg);

  double mme_err(0.0);
#endif

  // explicit time integration
  int nstep(0), nstep_dbg(0);
  double gcl, gcl_err(0.0);
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

      if (dim == 3) {
        cn *= 4.0 / 6.0;

        WhetStone::VectorPolynomial cn1, cn2;
        maps->NansonFormula(f, t, vec_vel[f], cn1);
        maps->NansonFormula(f, t + dt, vec_vel[f], cn2);

        cn += (1.0 / 6.0) * (cn1 + cn2);
      }

      (*vel)[f] = vec_vel[f] * cn;
    }

    // calculate cell velocity at time t+dt/2 and change its sign
    Entity_ID_List faces;
    std::vector<int> dirs;
    WhetStone::MatrixPolynomial C;
    WhetStone::VectorPolynomial tmp;

    for (int c = 0; c < ncells_owned; ++c) {
      mesh0->cell_get_faces(c, &faces);

      std::vector<WhetStone::VectorPolynomial> vvf;
      for (int n = 0; n < faces.size(); ++n) {
        vvf.push_back(vec_vel[faces[n]]);
      }

      maps->VelocityCell(c, vvf, tmp);
      maps->Cofactors(c, t + dt/2, tmp, C);
      (*cell_vel)[c].Multiply(C, tmp, true);

      // selecting new mean value
      std::vector<const WhetStone::Polynomial*> nvf;
      if (maps_name == "VEM") {
        for (int n = 0; n < faces.size(); ++n) {
          nvf.push_back(&(*vel)[faces[n]]);
        }
        dg.CoVelocityCell(c, nvf, (*cell_vel)[c]);
      }

      for (int i = 0; i < dim; ++i) {
        (*cell_vel)[c][i] *= -1.0;
      }

#ifdef DEBUG
      maps->VelocityCell(c, vvf, tmp);
      auto maps_dbg = std::make_shared<WhetStone::MeshMaps_FEM>(mesh0, mesh1);
      maps_dbg->Cofactors(c, t + dt/2, tmp, C);
      (*cell_vel_dbg)[c].Multiply(C, tmp, true);
      dg.CoVelocityCell(c, nvf, (*cell_vel_dbg)[c]);
#endif
    }

    // calculate determinant of Jacobian at time t+dt
    for (int c = 0; c < ncells_owned; ++c) {
      mesh0->cell_get_faces(c, &faces);
      std::vector<WhetStone::VectorPolynomial> vvf;

      for (int n = 0; n < faces.size(); ++n) {
        vvf.push_back(vec_vel[faces[n]]);
      }

      maps->JacobianDet(c, t, vvf, (*jac0)[c]);
      maps->JacobianDet(c, t + dt, vvf, (*jac1)[c]);

#ifdef DEBUG
      // maps->VelocityCell(c, vvf, tmp);
      // std::shared_ptr<WhetStone::MeshMaps_VEM> maps_aux = std::make_shared<WhetStone::MeshMaps_VEM>(mesh0, mesh1);
      // maps_aux->JacobianDet(t, tmp, (*jac0)[c]);
      // maps_aux->JacobianDet(t + dt, tmp, (*jac1)[c]);
      std::shared_ptr<WhetStone::MeshMaps> maps_dbg = std::make_shared<WhetStone::MeshMaps_FEM>(mesh0, mesh1);
      maps_dbg->JacobianDet(c, t, vvf, (*jac0_dbg)[c]);
#endif
    }

    // statistics: GCL
    for (int c = 0; c < ncells_owned; ++c) {
      double vol = mesh0->cell_volume(c);
      mesh0->cell_get_faces_and_dirs(c, &faces, &dirs);
      gcl = ((*jac1)[c](0, 0) - (*jac0)[c](0, 0)) * vol;

      for (int n = 0; n < faces.size(); ++n) {
        int f = faces[n];
        gcl -= (*vel)[f].Value(mesh0->face_centroid(f)) * dirs[n] * dt;
      }
      gcl_err += (gcl * gcl) * vol; 
    }

    // populate operators
    op->UpdateMatrices(vel.ptr());
    op_adv->UpdateMatrices();
    op_reac0->UpdateMatrices(p1.ptr());
    op_reac1->UpdateMatrices(p1.ptr());

#ifdef DEBUG
double dtref(1.0 / 10);
int it = (t + 1e-10) / dtref;
if(fabs(it * dtref - t) < 1e-12) {
    double dot0, dot1;
    CompositeVector dbg(*global_reac0->rhs());
    global_reac0->Apply(*p1_dbg, dbg);
    p1_dbg->Dot(dbg, &dot0);

    op_reac0_dbg->UpdateMatrices(p1_dbg.ptr());
    op_reac0_dbg->global_operator()->Apply(*p1_dbg, dbg);
    p1_dbg->Dot(dbg, &dot1);

    mme_err += fabs(dot0 - dot1) / dot0;

    dbg.PutScalar(0.0);
    op_adv_dbg0->UpdateMatrices();
    op_adv_dbg0->global_operator()->Apply(*p1_dbg, dbg);
    p1_dbg->Dot(dbg, &dot0);

    op_adv_dbg1->UpdateMatrices();
    op_adv_dbg1->global_operator()->Apply(*p1_dbg, dbg);
    p1_dbg->Dot(dbg, &dot1);

    ame_err += fabs(dot0 + dot1);
    nstep_dbg++;
}
#endif

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
    nstep++;
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

  err_tmp = gcl_err / nstep;
  mesh1->get_comm()->SumAll(&err_tmp, &gcl_err, 1);

  // error tests
  pl2_err = std::pow(pl2_err, 0.5);
  CHECK(pl2_err < 0.12 / (order + 1));

  gcl_err = std::pow(gcl_err, 0.5);

  if (MyPID == 0) {
    printf("L2(p0)=%12.8g  Inf(p0)=%12.8g  dMass=%12.8g  GCL=%12.8g  Err(area)=%12.8g\n", 
        pl2_err, pinf_err, mass1 - mass0, gcl_err, 1.0 - area);
  }

#ifdef DEBUG
  if (MyPID == 0)
    printf("L2 jacobian quadrature=%12.8g, advection quadrature=%12.8g\n", mme_err / nstep_dbg, ame_err / nstep_dbg);
#endif

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


const int N = 1;

TEST(REMAP2D_DG0_DUAL_VEM) {
  // RemapTestsDual(2, 0, "dg modal", "VEM", 16 * N, 16 * N, 0.05 / N);
  RemapTestsDual(2, 0, "dg modal", "VEM", 0, 0, 0.05 / N);
}

TEST(REMAP2D_DG1_DUAL_VEM) {
  // RemapTestsDual(2, 1, "dg modal", "VEM", 16 * N, 16 * N, 0.05 / N);
  RemapTestsDual(2, 1, "dg modal", "VEM", 0, 0, 0.05 / N);
}

TEST(REMAP3D_DG0_DUAL_VEM) {
  RemapTestsDual(3, 0, "dg modal", "VEM", 5, 5, 0.2);
}


/*
TEST(REMAP3D_DG1_DUAL_VEM) {
  RemapTestsDual(3, 1, "dg modal", "VEM", 5, 5, 0.1);
}
*/


/*
TEST(REMAP2D_DG1_QUADRATURE_ERROR) {
  RemapTestsDual(2, 1, "dg modal", "VEM", 16, 16, 0.1);
  RemapTestsDual(2, 1, "dg modal", "VEM", 16 *  2, 16 *  2, 0.1 / 2);
  RemapTestsDual(2, 1, "dg modal", "VEM", 16 *  4, 16 *  4, 0.1 / 4);
  RemapTestsDual(2, 1, "dg modal", "VEM", 16 *  8, 16 *  8, 0.1 / 8);
  RemapTestsDual(2, 1, "dg modal", "VEM", 16 * 16, 16 * 16, 0.1 / 16);
  RemapTestsDual(2, 1, "dg modal", "VEM", 16 * 32, 16 * 32, 0.1 / 32);
}
*/
