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
#include "LinearOperatorPCG.hh"
#include "MeshFactory.hh"
#include "MFD3D_Electromagnetics.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Accumulation.hh"
#include "ElectromagneticsMHD.hh"
#include "ElectromagneticsMHD_TM.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

#include "AnalyticElectromagnetics04.hh"
#include "AnalyticElectromagnetics05.hh"

/* *****************************************************************
* Testing operators for Maxwell-type problems: 2D 
* Magnetic flux B = (Bx, By, 0), electrif field E = (0, 0, Ez)
***************************************************************** */
template<class Analytic>
void ResistiveMHD2D(double dt, double tend, 
                    int nx, int ny, 
                    double Xa, double Ya, double Xb, double Yb,
                    const std::string& name) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Resistive MHD, TM mode, dt=" 
                            << dt << ", name: " << name << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_electromagnetics_TM.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a MSTK mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK}));

  RCP<const Mesh> mesh;
  if (name == "structured") {
    mesh = meshfactory(Xa, Ya, Xb, Yb, nx, ny, gm, true, true);
  } else {
    mesh = meshfactory(name, gm, true, true);
  }

  // create resistivity coefficient
  double told(0.0), tnew(dt);
  Analytic ana(mesh);
  WhetStone::Tensor Kc(2, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, tnew);
    K->push_back(Kc);
  }

  // create miscalleneous data
  int nnodes_owned = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Teuchos::RCP<BCs> bc1 = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, SCHEMA_DOFS_SCALAR));

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("electromagnetics operator");
  Teuchos::RCP<ElectromagneticsMHD_TM> op_mhd = Teuchos::rcp(new ElectromagneticsMHD_TM(olist, mesh));
  op_mhd->SetBCs(bc1, bc1);

  // create/extract solution maps
  Teuchos::RCP<Operator> global_op = op_mhd->global_operator();
  const CompositeVectorSpace& cvs_e = global_op->DomainMap();

  Teuchos::RCP<CompositeVectorSpace> cvs_b = Teuchos::rcp(new CompositeVectorSpace());
  cvs_b->SetMesh(mesh)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector E(cvs_e);
  CompositeVector B(*cvs_b);

  // set up initial guess for a time-dependent problem
  Epetra_MultiVector& Ee = *E.ViewComponent("node");
  Epetra_MultiVector& Bf = *B.ViewComponent("face");

  AmanziGeometry::Point xv(2);
  std::vector<int> edirs;
  AmanziMesh::Entity_ID_List cells, faces;

  Ee.PutScalar(0.0);
  Bf.PutScalar(0.0);

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh->node_get_coordinates(v, &xv);
    Ee[0][v] = (ana.electric_exact(xv, told))[2];
  }

  for (int f = 0; f < nfaces_owned; ++f) {
    double area = mesh->face_area(f);
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);

    Bf[0][f] = (normal * ana.magnetic_exact(xf, told)) / area;
  } 

  // CompositeVector B0(B);
  // Epetra_MultiVector& B0f = *B0.ViewComponent("face");

  int cycle(0);
  double energy0(1e+99), divB0(0.0);
  while (told + dt/2 < tend) {
    // set up the diffusion operator
    global_op->Init();
    op_mhd->SetTensorCoefficient(K);
    op_mhd->UpdateMatrices();

    // Add an accumulation term using dt=1 since time step is taken into
    // account in the system modification routine. Kc=constant FIXME
    CompositeVector phi(cvs_e);
    phi.PutScalar(1.0 / Kc(0,0));

    Teuchos::RCP<Accumulation> op_acc = Teuchos::rcp(new Accumulation(AmanziMesh::NODE, global_op));
    op_acc->SetBCs(bc1, bc1);
    op_acc->AddAccumulationTerm(phi, 1.0, "node");

    // update BCs
    std::vector<int>& bc_model = bc1->bc_model();
    std::vector<double>& bc_value = bc1->bc_value();

    for (int v = 0; v < nnodes_wghost; ++v) {
      mesh->node_get_coordinates(v, &xv);

      if (fabs(xv[0] - Xa) < 1e-6 || fabs(xv[0] - Xb) < 1e-6 ||
          fabs(xv[1] - Ya) < 1e-6 || fabs(xv[1] - Yb) < 1e-6) {
        bc_model[v] = OPERATOR_BC_DIRICHLET;
        bc_value[v] = (ana.electric_exact(xv, tnew - dt/2))[2];
      }
    }

    // apply BCs and assemble
    op_mhd->ModifyMatrices(E, B, dt);
    op_mhd->ApplyBCs(true, true);
    op_acc->ApplyBCs();
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();

    ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
    global_op->InitPreconditioner("Hypre AMG", slist);

    // Solve the problem.
    ParameterList lop_list = plist.sublist("solvers").sublist("silent").sublist("pcg parameters");
    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        solver(global_op, global_op);
    solver.Init(lop_list);

    CompositeVector& rhs = *global_op->rhs();
    int ierr = solver.ApplyInverse(rhs, E);

    double energy = op_mhd->CalculateMagneticEnergy(B);
    op_mhd->ModifyFields(E, B, dt);

    // CHECK(energy < energy0);
    energy0 = energy;

    cycle++;
    told = tnew;
    tnew += dt;

    // reconstruction
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 2);

    CompositeVector Bvec(*cvs);
    Epetra_MultiVector& sol = *Bvec.ViewComponent("cell"); 
    sol.PutScalar(0.0);

    std::vector<int> dirs;
    AmanziMesh::Entity_ID_List faces;

    double avgB(0.0), divB(0.0), errB(0.0);
    for (int c = 0; c < ncells_owned; ++c) {
      double vol = mesh->cell_volume(c);
      const Amanzi::AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      double tmp(0.0);
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        double area = mesh->face_area(f);
        const Amanzi::AmanziGeometry::Point& xf = mesh->face_centroid(f);
        for (int k = 0; k < 2; ++k) {
          sol[k][c] += Bf[0][f] * dirs[n] * area * (xf[k] - xc[k]) / vol;
        }        
        tmp += Bf[0][f] * dirs[n] * area / vol;
      }
      avgB += std::fabs(sol[0][c]);
      divB += tmp * tmp * vol; 
      errB += vol * (std::pow(sol[0][c], 2.0) + std::pow(sol[1][c], 2.0));
    }

    if (cycle == 1) divB0 = divB;
    CHECK_CLOSE(divB0, divB, 1e-8);

    if (MyPID == 0) {
      std::cout << "time: " << told << "  ||r||=" << solver.residual() 
                << " itr=" << solver.num_itrs() << "  energy= " << energy 
                << "  avgB=" << avgB / ncells_owned 
                << "  divB=" << std::pow(divB, 0.5) 
                << "  ||B||=" << std::pow(errB, 0.5) << std::endl;
    }

    // visualization
    if (MyPID == 0) {
      GMV::open_data_file(*mesh, "operators.gmv");
      GMV::start_data();
      GMV::write_cell_data(sol, 0, "Bx");
      GMV::write_cell_data(sol, 1, "By");
      GMV::write_node_data(Ee, 0, "Ez");
      GMV::close_data_file();
    }
  }

  // compute electric and magnetic errors
  double enorm, el2_err, einf_err;
  double bnorm, bl2_err, binf_err;
  ana.ComputeNodeError(Ee, told - dt/2, enorm, el2_err, einf_err);
  ana.ComputeFaceError(Bf, told, bnorm, bl2_err, binf_err);

  if (MyPID == 0) {
    if (enorm != 0.0) el2_err /= enorm; 
    if (bnorm != 0.0) bl2_err /= bnorm; 
    printf("L2(e)=%10.7f  Inf(e)=%9.6f  L2(b)=%10.7f  Inf(b)=%9.6f\n",
        el2_err, einf_err, bl2_err, binf_err);
    // CHECK(el2_err < tolerance);
  }
}


/* *****************************************************************
* Testing operators for Maxwell-type problems: 3D
* **************************************************************** */
template<class Analytic>
void ResistiveMHD3D(double dt, double tend, bool convergence,
                    int nx, int ny, int nz,
                    double Xa, double Ya, double Za, double Xb, double Yb, double Zb,
                    const std::string& name) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Resistive MHD 3D, dt=" 
                            << dt << ", name: " << name << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_electromagnetics.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a MSTK mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, &comm));

  MeshFactory meshfactory(&comm);
  meshfactory.preference(FrameworkPreference({MSTK}));

  bool request_faces(true), request_edges(true);
  RCP<const Mesh> mesh;
  if (name == "structured") {
    mesh = meshfactory(Xa, Ya, Za, Xb, Yb, Zb, nx, ny, nz, gm, request_faces, request_edges);
  } else {
    mesh = meshfactory(name, gm, request_faces, request_edges);
    // mesh = meshfactory("test/hex_split_faces5.exo", gm, request_faces, request_edges);
  }

  // create resistivity coefficient
  double told(0.0), tnew(dt);
  Analytic ana(mesh);
  WhetStone::Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, tnew);
    K->push_back(Kc);
  }

  // create boundary data
  int nedges_owned = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::OWNED);
  int nedges_wghost = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::USED);

  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Teuchos::RCP<BCs> bc1 = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, SCHEMA_DOFS_SCALAR));
  Teuchos::RCP<BCs> bc2 = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, SCHEMA_DOFS_SCALAR));

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("electromagnetics operator");
  Teuchos::RCP<ElectromagneticsMHD> op_mhd = Teuchos::rcp(new ElectromagneticsMHD(olist, mesh));
  op_mhd->SetBCs(bc1, bc1);
  if (!convergence) op_mhd->AddBCs(bc2, bc2);

  // create/extract solution maps
  Teuchos::RCP<Operator> global_op = op_mhd->global_operator();
  const CompositeVectorSpace& cvs_e = global_op->DomainMap();

  Teuchos::RCP<CompositeVectorSpace> cvs_b = Teuchos::rcp(new CompositeVectorSpace());
  cvs_b->SetMesh(mesh)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector E(cvs_e);
  CompositeVector B(*cvs_b);

  // set up initial guess for a time-dependent problem
  Epetra_MultiVector& Ee = *E.ViewComponent("edge");
  Epetra_MultiVector& Bf = *B.ViewComponent("face");

  Ee.PutScalar(0.0);
  Bf.PutScalar(0.0);

  for (int e = 0; e < nedges_owned; ++e) {
    double len = mesh->edge_length(e);
    const AmanziGeometry::Point& tau = mesh->edge_vector(e);
    const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

    Ee[0][e] = (ana.electric_exact(xe, told) * tau) / len;
  }

  for (int f = 0; f < nfaces_owned; ++f) {
    double area = mesh->face_area(f);
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);

    Bf[0][f] = (ana.magnetic_exact(xf, told) * normal) / area;
  }
  // CompositeVector B0(B);
  // Epetra_MultiVector& B0f = *B0.ViewComponent("face");

  int cycle(0);
  double energy0(1e+99), divB0(0.0);
  while (told + dt/2 < tend) {
    // set up the diffusion operator
    global_op->Init();
    op_mhd->SetTensorCoefficient(K);
    op_mhd->UpdateMatrices();

    // update BCs
    std::vector<int>& bc_model = bc1->bc_model();
    std::vector<double>& bc_value = bc1->bc_value();

    std::vector<int>& bc_model2 = bc2->bc_model();
    std::vector<Point>& bc_value2 = bc2->bc_value_point();

    std::vector<int> edirs;
    AmanziMesh::Entity_ID_List cells, edges;

    for (int f = 0; f < nfaces_wghost; ++f) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);

      if (fabs(xf[0] - Xa) < 1e-6 && !convergence) {
        bc_model2[f] = OPERATOR_BC_NEUMANN;
        bc_value2[f] = Point(0.0, 0.0, 1.0);
      } 
      else if ((fabs(xf[0] - Xa) < 1e-6 && convergence) || fabs(xf[0] - Xb) < 1e-6 ||
                fabs(xf[1] - Ya) < 1e-6 || fabs(xf[1] - Yb) < 1e-6 ||
                fabs(xf[2] - Za) < 1e-6 || fabs(xf[2] - Zb) < 1e-6) {
        mesh->face_get_edges_and_dirs(f, &edges, &edirs);
        int nedges = edges.size();
        for (int i = 0; i < nedges; ++i) {
          int e = edges[i];
          double len = mesh->edge_length(e);
          const AmanziGeometry::Point& tau = mesh->edge_vector(e);
          const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

          bc_model[e] = OPERATOR_BC_DIRICHLET;
          bc_value[e] = (ana.electric_exact(xe, tnew - dt/2) * tau) / len;
        }
      }
    }

    // BCs, sources, and assemble
    op_mhd->ModifyMatrices(E, B, dt);
    op_mhd->ApplyBCs(true, true);
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();

    ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
    global_op->InitPreconditioner("Hypre AMG", slist);

    // Solve the problem.
    ParameterList lop_list = plist.sublist("solvers").sublist("silent").sublist("pcg parameters");
    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        solver(global_op, global_op);
    solver.Init(lop_list);

    CompositeVector& rhs = *global_op->rhs();
    int ierr = solver.ApplyInverse(rhs, E);

    double energy = op_mhd->CalculateMagneticEnergy(B);
    op_mhd->ModifyFields(E, B, dt);

    if (!convergence) CHECK(energy < energy0);
    energy0 = energy;

    cycle++;
    told = tnew;
    tnew += dt;

    // reconstruction
    // -- magnetic field
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 3);

    CompositeVector Bvec(*cvs);
    Epetra_MultiVector& sol_b = *Bvec.ViewComponent("cell"); 
    sol_b.PutScalar(0.0);

    std::vector<int> dirs;
    AmanziMesh::Entity_ID_List faces;

    double avgB(0.0), divB(0.0), errB(0.0);
    for (int c = 0; c < ncells_owned; ++c) {
      double vol = mesh->cell_volume(c);
      const Amanzi::AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      double tmp(0.0);
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        double area = mesh->face_area(f);
        const Amanzi::AmanziGeometry::Point& xf = mesh->face_centroid(f);
        for (int k = 0; k < 3; ++k) {
          sol_b[k][c] += Bf[0][f] * dirs[n] * area * (xf[k] - xc[k]) / vol;
        }        
        tmp += Bf[0][f] * dirs[n] * area / vol;
      }
      avgB += std::fabs(sol_b[0][c]);
      divB += tmp * tmp * vol; 

      AmanziGeometry::Point v1(ana.magnetic_exact(xc, told));
      AmanziGeometry::Point v2(sol_b[0][c], sol_b[1][c], sol_b[2][c]);
      errB += vol * (L22(v1 - v2));
    }

    // -- electric field
    CompositeVector Evec(*cvs);
    Epetra_MultiVector& sol_e = *Evec.ViewComponent("cell"); 
    sol_e.PutScalar(0.0);

    double avgE(0.0);
    WhetStone::MFD3D_Electromagnetics mfd(mesh); 
    WhetStone::Tensor Ic(3, 1);
    Ic(0, 0) = 1.0;

    for (int c = 0; c < ncells_owned; ++c) {
      mesh->cell_get_edges(c, &edges);
      int nedges = edges.size();

      WhetStone::DenseMatrix R(nedges, 3), W(nedges, nedges);
      WhetStone::DenseVector v1(nedges), v2(3);

      for (int n = 0; n < nedges; ++n) {
        v1(n) = Ee[0][edges[n]];
      }

      mfd.L2consistencyInverse(c, Ic, R, W, true);
      R.Multiply(v1, v2, true);

      double vol = mesh->cell_volume(c);
      for (int k = 0; k < 3; ++k) {
        sol_e[k][c] = v2(k) / vol;
      }
    }


    if (cycle == 1) divB0 = divB;
    CHECK_CLOSE(divB0, divB, 1e-8);

    if (MyPID == 0) {
      std::cout << "time: " << told << "  ||r||=" << solver.residual() 
                << " itr=" << solver.num_itrs() << "  energy= " << energy 
                << "  avgB=" << avgB / ncells_owned 
                << "  divB=" << std::pow(divB, 0.5) 
                << "  errB=" << std::pow(errB, 0.5) << std::endl;
    }

    // visualization
    if (MyPID == 0 && (cycle % 5 == 0)) {
      GMV::open_data_file(*mesh, "operators.gmv");
      GMV::start_data();
      GMV::write_cell_data(sol_b, 0, "Bx");
      GMV::write_cell_data(sol_b, 1, "By");
      GMV::write_cell_data(sol_b, 2, "Bz");

      GMV::write_cell_data(sol_e, 0, "Ex");
      GMV::write_cell_data(sol_e, 1, "Ey");
      GMV::write_cell_data(sol_e, 2, "Ez");
      GMV::close_data_file();
    }
  }

  // compute electric and magnetic errors
  if (convergence) {
    double enorm, el2_err, einf_err;
    double bnorm, bl2_err, binf_err;
    ana.ComputeEdgeError(Ee, told - dt/2, enorm, el2_err, einf_err);
    ana.ComputeFaceError(Bf, told, bnorm, bl2_err, binf_err);

    if (MyPID == 0) {
      if (enorm != 0.0) el2_err /= enorm; 
      if (bnorm != 0.0) bl2_err /= bnorm; 
      printf("L2(e)=%10.7f  Inf(e)=%9.6f  L2(b)=%10.7f  Inf(b)=%9.6f\n",
          el2_err, einf_err, bl2_err, binf_err);
    }
  }
}


TEST(RESISTIVE_MHD2D_RELAX) {
  ResistiveMHD2D<AnalyticElectromagnetics04>(0.7, 5.9, 8,18, -4.0,-10.0, 4.0,10.0, "structured");
}

TEST(RESISTIVE_MHD2D_CONVERGENCE) {
  ResistiveMHD2D<AnalyticElectromagnetics05>(0.01, 0.1, 10,10, 0.0,0.0, 1.0,1.0, "test/random10.exo");
  // ResistiveMHD2D<AnalyticElectromagnetics05>(0.01 / 2, 0.1, 20,20, 0.0,0.0, 1.0,1.0, "structured");
}

TEST(RESISTIVE_MHD3D_RELAX) {
  ResistiveMHD3D<AnalyticElectromagnetics04>(0.1, 0.5, false, 10,8,22, -4.0,-4.0,-10.0, 4.0,4.0,10.0, "structured");
}

TEST(RESISTIVE_MHD3D_CONVERGENCE) {
  // ResistiveMHD3D<AnalyticElectromagnetics05>(0.01, 0.1, true, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "structured");
  ResistiveMHD3D<AnalyticElectromagnetics05>(0.01, 0.1, true, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "test/kershaw08.exo");
  // ResistiveMHD3D<AnalyticElectromagnetics05>(0.01 / 2, 0.1, true, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "test/random3D_20.exo");
}

