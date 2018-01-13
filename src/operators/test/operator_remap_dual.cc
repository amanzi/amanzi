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
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "WhetStone_typedefs.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

#include "AnalyticDG00.hh"
#include "AnalyticDG04.hh"
#include "AnalyticDG05.hh"


/* *****************************************************************
* Remap of polynomilas in two dimensions. Explicit time scheme.
* Dual formulation places gradient and jumps on a test function.
***************************************************************** */
void RemapTestsDual(int dim, int order_p, int order_u,
                    std::string maps_name,
                    int nx, int ny, double dt) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: " << dim << "D remap, dual formulation:"
                            << " mesh=" << ((nx == 0) ? "polygonal" : "square")
                            << ", orders: " << order_p << " " << order_u 
                            << ", maps=" << maps_name << std::endl;

  // polynomial space
  WhetStone::Polynomial pp(dim, order_p);
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
  Entity_ID_List nodeids, faces;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh1->node_get_coordinates(v, &xv);
    yv = xv;

    double ds(0.0001);
    for (int i = 0; i < 10000; ++i) {
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
    std::shared_ptr<WhetStone::MeshMaps_VEM> maps_vem;
    maps_vem = std::make_shared<WhetStone::MeshMaps_VEM>(mesh0, mesh1);
    maps_vem->set_order(order_u);  // higher-order velocity reconstruction
    maps = maps_vem;
  }

  // numerical integration
  WhetStone::NumericalIntegration numi(mesh0);

  // create and initialize cell-based field 
  CompositeVectorSpace cvs1;
  cvs1.SetMesh(mesh0)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, nk);
  Teuchos::RCP<CompositeVector> p1 = Teuchos::rcp(new CompositeVector(cvs1));
  Epetra_MultiVector& p1c = *p1->ViewComponent("cell", true);

  // we need dg to compute scaling of basis functions
  WhetStone::DG_Modal dg(order_p, mesh0);

  AnalyticDG04 ana(mesh0, order_p);

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh0->cell_centroid(c);
    WhetStone::Polynomial coefs;
    ana.TaylorCoefficients(xc, 0.0, coefs);

    WhetStone::DenseVector data;
    coefs.GetPolynomialCoefficients(data);

    double a, b;
    for (auto it = coefs.begin(); it.end() <= coefs.end(); ++it) {
      int n = it.PolynomialPosition();

      dg.TaylorBasis(c, it, &a, &b);
      p1c[n][c] = data(n) / a;
      p1c[0][c] += data(n) * b;
    }
  }

  // initial mass
  double mass0(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    std::vector<double> data;
    for (int i = 0; i < nk; ++i) {
      data.push_back(p1c[i][c]);
    }
    WhetStone::Polynomial poly(dg.CalculatePolynomial(c, data));
    mass0 += numi.IntegratePolynomialCell(c, poly);
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
  plist.set<std::string>("method", "dg modal")
       .set<int>("method order", order_p)
       .set<bool>("jump operator on test function", true);

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
       .set<bool>("gradient operator on test function", true);
  plist.sublist("schema domain").set<std::string>("base", "cell");
  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<PDE_Abstract> op_adv = Teuchos::rcp(new PDE_Abstract(plist, global_op));

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

  // explicit time integration
  int nstep(0), nstep_dbg(0);
  double gcl, gcl_err(0.0);
  double t(0.0), tend(1.0);
  while(t < tend - dt/2) {
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

    // calculate various geometric quantities in reference framework
    WhetStone::MatrixPolynomial C;

    for (int c = 0; c < ncells_owned; ++c) {
      // product -C^t u at time t+dt/2
      maps->Cofactors(t + dt/2, J[c], C);
      (*cell_vel)[c].Multiply(C, uc[c], true);

      for (int i = 0; i < dim; ++i) {
        (*cell_vel)[c][i] *= -1.0;
      }

      // determinant of Jacobian
      maps->Determinant(t, J[c], (*jac0)[c]);
      maps->Determinant(t + dt, J[c], (*jac1)[c]);
    }

    // statistics: GCL
    std::vector<int> dirs;

    for (int c = 0; c < ncells_owned; ++c) {
      double vol = mesh0->cell_volume(c);
      mesh0->cell_get_faces_and_dirs(c, &faces, &dirs);
      gcl = numi.IntegratePolynomialCell(c, (*jac1)[c]) 
          - numi.IntegratePolynomialCell(c, (*jac0)[c]);

      for (int n = 0; n < faces.size(); ++n) {
        int f = faces[n];
        gcl -= (*vel)[f].Value(mesh0->face_centroid(f)) * dirs[n] * dt;
      }
      gcl_err = std::max(gcl_err, fabs(gcl) / mesh0->cell_volume(c)); 
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
    criteria.push_back("relative rhs");
    plist.set<double>("error tolerance", 1e-12)
         .set<Teuchos::Array<std::string> >("convergence criteria", criteria)
         .sublist("verbose object")
         .set<std::string>("verbosity level", "none");
    pcg.Init(plist);
    pcg.ApplyInverse(g, p2);

    if (MyPID == 0)
      printf("time=%8.5f  pcg=(%d, %9.4g)\n", t, pcg.num_itrs(), pcg.residual());

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
  double pl2_err(0.0), pinf_err(0.0), area(0.0);
  double mass1(0.0), ql2_err(0.0), qinf_err(0.0);

  Entity_ID_List nodes;
  std::vector<int> dirs;
  AmanziGeometry::Point v0(dim), v1(dim), tau(dim);

  for (int c = 0; c < ncells_owned; ++c) {
    double area_c(mesh1->cell_volume(c));

    std::vector<double> data;
    for (int i = 0; i < nk; ++i) {
      data.push_back(p2c[i][c]);
    }
    WhetStone::Polynomial poly(dg.CalculatePolynomial(c, data));

    if (nk == 1) {
      // const AmanziGeometry::Point& xg = maps->cell_geometric_center(1, c);
      const AmanziGeometry::Point& xg = mesh1->cell_centroid(c);
      double tmp = p2c[0][c] - ana.function_exact(xg, 0.0);

      pinf_err = std::max(pinf_err, fabs(tmp));
      pl2_err += tmp * tmp * area_c;
    }
    else {
      mesh0->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();  
      for (int i = 0; i < nnodes; ++i) {
        mesh0->node_get_coordinates(nodes[i], &v0);
        mesh1->node_get_coordinates(nodes[i], &v1);

        double tmp = poly.Value(v0);
        tmp -= ana.function_exact(v1, 0.0);
        pinf_err = std::max(pinf_err, fabs(tmp));
        pl2_err += tmp * tmp * area_c / nnodes;
      }
    }

    area += area_c;

    WhetStone::Polynomial tmp((*jac1)[c]);
    tmp.ChangeOrigin(mesh0->cell_centroid(c));
    poly *= tmp;

    double mass1_c = numi.IntegratePolynomialCell(c, poly);
    mass1 += mass1_c;

    // optional projection on the space of polynomials 
    if (order_p > 0 && dim == 2) {
      poly = dg.CalculatePolynomial(c, data);

      mesh0->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();  

      WhetStone::VectorPolynomial vc(dim, 0);
      std::vector<WhetStone::VectorPolynomial> vvf;

      for (int i = 0; i < nfaces; ++i) {
        int f = faces[i];
        const AmanziGeometry::Point& xf = mesh1->face_centroid(f);
        mesh0->face_get_nodes(f, &nodes);

        mesh0->node_get_coordinates(nodes[0], &v0);
        mesh0->node_get_coordinates(nodes[1], &v1);
        double f0 = poly.Value(v0);
        double f1 = poly.Value(v1);

        WhetStone::VectorPolynomial vf(dim - 1, 1);
        vf[0].Reshape(dim - 1, order_p);
        if (order_p == 1) {
          vf[0](0, 0) = (f0 + f1) / 2; 
          vf[0](1, 0) = f1 - f0; 
        } else if (order_p == 2) {
          double f2 = poly.Value((v0 + v1) / 2);
          vf[0](0, 0) = f2;
          vf[0](1, 0) = f1 - f0;
          vf[0](2, 0) = -4 * f2 + 2 * f0 + 2 * f1;
        } else {
          ASSERT(0);
        }

        mesh1->node_get_coordinates(nodes[0], &v0);
        mesh1->node_get_coordinates(nodes[1], &v1);

        std::vector<AmanziGeometry::Point> tau;
        tau.push_back(v1 - v0);
        vf[0].InverseChangeCoordinates(xf, tau);
        vvf.push_back(vf);
      }

      auto moments = std::make_shared<WhetStone::DenseVector>();
      if (order_p == 2) {
        moments->Reshape(1);
        (*moments)(0) = mass1_c / mesh1->cell_volume(c);
      }

      WhetStone::ProjectorH1 projector(mesh1);
      projector.EllipticCell_Pk(c, order_p, vvf, moments, vc);
      vc[0].ChangeOrigin(mesh1->cell_centroid(c));

      if (order_p == 1) {
        vc[0](0, 0) = mass1_c / mesh1->cell_volume(c);
      }

      // error in the projected solution vc[0]
      mesh0->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();  
      for (int i = 0; i < nnodes; ++i) {
        mesh1->node_get_coordinates(nodes[i], &v1);

        double tmp = vc[0].Value(v1);
        tmp -= ana.function_exact(v1, 0.0);
        qinf_err = std::max(qinf_err, fabs(tmp));
        ql2_err += tmp * tmp * area_c / nnodes;
      }
    }
  }

  // parallel collective operations
  double err_in[4] = {pl2_err, area, mass1, ql2_err};
  double err_out[4];
  mesh1->get_comm()->SumAll(err_in, err_out, 4);

  double err_tmp = pinf_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &pinf_err, 1);

  err_tmp = qinf_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &qinf_err, 1);

  err_tmp = gcl_err;
  mesh1->get_comm()->MaxAll(&err_tmp, &gcl_err, 1);

  // error tests
  pl2_err = std::pow(err_out[0], 0.5);
  ql2_err = std::pow(err_out[3], 0.5);
  CHECK(pl2_err < 0.12 / (order_p + 1));

  if (MyPID == 0) {
    printf("nx=%3d  L2=%12.8g %12.8g  Inf=%12.8g %12.8f dMass=%10.4g  GCL_Inf=%10.6g  dArea=%10.6g\n", 
        nx, pl2_err, ql2_err, pinf_err, qinf_err, err_out[2] - mass0, gcl_err, 1.0 - err_out[1]);
  }

  // visualization
  if (MyPID == 0) {
    GMV::open_data_file(*mesh1, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p2c, 0, "remaped");
    if (order_p > 0) {
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
const double q = 2.0;  // 2.82842712474619;

// TEST(REMAP2D_DG1_DUAL_VEM_TMP) {
//   RemapTestsDual(2, 1, 1, "VEM", 0, 0, 0.05 / N);
//   RemapTestsDual(2, 1, 2, "VEM", 0, 0, 0.05 / N);
// }

TEST(REMAP2D_DG0_DUAL_VEM) {
  // RemapTestsDual(2, 0, 1, "VEM", 16 * N, 16 * N, 0.05 / N);
  RemapTestsDual(2, 0, 1, "VEM", 0, 0, 0.05 / N);
}

TEST(REMAP2D_DG1_DUAL_VEM) {
  // RemapTestsDual(2, 1, 2, "VEM", 16 * N, 16 * N, 0.05 / N);
  RemapTestsDual(2, 1, 2, "VEM", 0, 0, 0.05 / N);
}

TEST(REMAP2D_DG2_DUAL_VEM) {
  // RemapTestsDual(2, 2, 3, "VEM", 10, 10, 0.05 / N);
  RemapTestsDual(2, 2, 3, "VEM", 0, 0, 0.05 / N);
}

TEST(REMAP3D_DG0_DUAL_VEM) {
  RemapTestsDual(3, 0, 1, "VEM", 5, 5, 0.2);
}

// TEST(REMAP3D_DG1_DUAL_VEM) {
//   RemapTestsDual(3, 1, 2, "VEM", 5, 5, 0.1);
// }

/*
TEST(REMAP2D_DG1_QUADRATURE_ERROR) {
  RemapTestsDual(2, 2, 2, "VEM", 16, 16, 0.05);
  RemapTestsDual(2, 2, 2, "VEM", 16 *  2, 16 *  2, 0.05 / q);
  RemapTestsDual(2, 2, 2, "VEM", 16 *  4, 16 *  4, 0.05 / std::pow(q, 2.0));
  RemapTestsDual(2, 2, 2, "VEM", 16 *  8, 16 *  8, 0.05 / std::pow(q, 3.0));
  RemapTestsDual(2, 2, 2, "VEM", 16 * 16, 16 * 16, 0.05 / std::pow(q, 4.0));
  // RemapTestsDual(2, 1, 1, "VEM", 16 * 32, 16 * 32, 0.05 / std::pow(q, 5.0));
}
*/
