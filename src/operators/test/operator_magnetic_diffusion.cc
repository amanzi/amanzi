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
#include "WhetStoneMeshUtils.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_MagneticDiffusion.hh"
#include "PDE_MagneticDiffusion_TM.hh"

#include "AnalyticElectromagnetics04.hh"
#include "AnalyticElectromagnetics05.hh"
#include "MeshDeformation.hh"

/* *****************************************************************
* Testing operators for Maxwell-type problems: 2D 
* Magnetic flux B = (Bx, By, 0), electric field E = (0, 0, Ez)
***************************************************************** */
template<class Analytic>
void MagneticDiffusion2D(double dt, double tend, 
                    int nx, int ny, 
                    double Xa, double Ya, double Xb, double Yb,
                    const std::string& name) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Magnetic diffusion, TM mode, dt=" 
                            << dt << ", name: " << name << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_electromagnetics_TM.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  RCP<const Mesh> mesh;
  if (name == "structured") {
    mesh = meshfactory.create(Xa, Ya, Xb, Yb, nx, ny, true, true);
  } else {
    mesh = meshfactory.create(name, true, true);
  }

  // create resistivity coefficient
  double told(0.0), tnew(dt);
  Analytic ana(mesh);
  WhetStone::Tensor Kc(2, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, tnew);
    K->push_back(Kc);
  }

  // create miscalleneous data
  int nnodes_owned = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nnodes_wghost = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<BCs> bc1 = Teuchos::rcp(new BCs(mesh, AmanziMesh::NODE, WhetStone::DOF_Type::SCALAR));

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("electromagnetics operator");
  Teuchos::RCP<PDE_MagneticDiffusion_TM> op_mag = Teuchos::rcp(new PDE_MagneticDiffusion_TM(olist, mesh));
  op_mag->SetBCs(bc1, bc1);

  // create/extract solution maps
  Teuchos::RCP<Operator> global_op = op_mag->global_operator();
  const CompositeVectorSpace& cvs_e = global_op->DomainMap();

  Teuchos::RCP<CompositeVectorSpace> cvs_b = Teuchos::rcp(new CompositeVectorSpace());
  cvs_b->SetMesh(mesh)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);

  CompositeVector E(cvs_e);
  CompositeVector B(*cvs_b);

  // set up initial guess for a time-dependent problem
  Epetra_MultiVector& Ee = *E.ViewComponent("node");
  Epetra_MultiVector& Bf = *B.ViewComponent("face");

  AmanziGeometry::Point xv(2);

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
    op_mag->SetTensorCoefficient(K);
    op_mag->UpdateMatrices();

    // Add an accumulation term using dt=1 since time step is taken into
    // account in the system modification routine. Kc=constant FIXME
    CompositeVector phi(cvs_e);
    phi.PutScalar(1.0 / Kc(0,0));

    Teuchos::RCP<PDE_Accumulation> op_acc = Teuchos::rcp(new PDE_Accumulation(AmanziMesh::NODE, global_op));
    op_acc->SetBCs(bc1, bc1);
    op_acc->AddAccumulationTerm(phi, 1.0, "node");

    // update electric boundary conditions
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
    op_mag->ModifyMatrices(E, B, dt);
    op_mag->ApplyBCs(true, true, true);
    op_acc->ApplyBCs();
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();

    ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
    global_op->InitializePreconditioner(slist);
    global_op->UpdatePreconditioner();

    // Solve the problem.
    ParameterList lop_list = plist.sublist("solvers").sublist("silent").sublist("pcg parameters");
    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        solver(global_op, global_op);
    solver.Init(lop_list);

    CompositeVector& rhs = *global_op->rhs();
    solver.ApplyInverse(rhs, E);

    double heat = op_mag->CalculateOhmicHeating(E);
    double energy = op_mag->CalculateMagneticEnergy(B);
    op_mag->ModifyFields(E, B, dt);

    CHECK(heat > 0.0);
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
    ana.GlobalOp("sum", &avgB, 1);
    ana.GlobalOp("sum", &divB, 1);
    ana.GlobalOp("sum", &errB, 1);

    if (cycle == 1) divB0 = divB;
    CHECK_CLOSE(divB0, divB, 1e-8);

    if (MyPID == 0) {
      std::cout << "time: " << told << "  ||r||=" << solver.residual() 
                << " itr=" << solver.num_itrs() << "  energy= " << energy 
                << "  heat= " << heat
                << "  avgB=" << avgB / Bf.GlobalLength()
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
    printf("nx=%d L2(e)=%10.7f  Inf(e)=%9.6f  L2(b)=%10.7f  Inf(b)=%9.6f\n",
        nx, el2_err, einf_err, bl2_err, binf_err);
    // CHECK(el2_err < tolerance);
  }
}


TEST(MAGNETIC_DIFFUSION2D_RELAX) {
  MagneticDiffusion2D<AnalyticElectromagnetics04>(0.7, 5.9, 8,18, -4.0,-10.0, 4.0,10.0, "structured");
}



/* *****************************************************************
* Testing operators for Maxwell-type problems: 3D
* **************************************************************** */
template<class Analytic>
void MagneticDiffusion3D(double dt, double tend, bool convergence,
                    int nx, int ny, int nz,
                    double Xa, double Ya, double Za, double Xb, double Yb, double Zb,
                    const std::string& name, int deform = 0) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Magnetic diffusion 3D, dt=" 
                            << dt << ", name: " << name << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_electromagnetics.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  bool request_faces(true), request_edges(true);
  RCP<Mesh> mesh;
  if (name == "structured")
    mesh = meshfactory.create(Xa, Ya, Za, Xb, Yb, Zb, nx, ny, nz, request_faces, request_edges);
  else
    mesh = meshfactory.create(name, request_faces, request_edges);
    // mesh = meshfactory.create("test/hex_split_faces5.exo", request_faces, request_edges);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  Analytic ana(mesh);

  if (deform > 0) {
    double vol0(0.0), vol1(0.0);
    for (int c = 0; c < ncells_owned; ++c) vol0 += mesh->cell_volume(c);
    DeformMesh(mesh, deform, 0.0);
    for (int c = 0; c < ncells_owned; ++c) vol1 += mesh->cell_volume(c);

    vol0 -= vol1;
    ana.GlobalOp("sum", &vol0, 1);
    if (MyPID == 0)
      std::cout << "volume change after deformation=" << vol0 << std::endl;
  }

  // create resistivity coefficient
  double told(0.0), tnew(dt);
  WhetStone::Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, tnew);
    K->push_back(Kc);
  }

  // create boundary data
  int nedges_owned = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Teuchos::RCP<BCs> bc1 = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, WhetStone::DOF_Type::SCALAR));
  Teuchos::RCP<BCs> bc2 = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("electromagnetics operator");
  Teuchos::RCP<PDE_MagneticDiffusion> op_mag = Teuchos::rcp(new PDE_MagneticDiffusion(olist, mesh));
  op_mag->SetBCs(bc1, bc1);
  if (!convergence) op_mag->AddBCs(bc2, bc2);

  // create/extract solution maps
  Teuchos::RCP<Operator> global_op = op_mag->global_operator();
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

  int cycle(0);
  double energy0(1e+99), divB0(0.0);
  while (told + dt/2 < tend) {
    // set up the diffusion operator
    global_op->Init();
    op_mag->SetTensorCoefficient(K);
    op_mag->UpdateMatrices();

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
    op_mag->ModifyMatrices(E, B, dt);
    op_mag->ApplyBCs(true, true, true);
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();

    ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
    global_op->InitializePreconditioner(slist);
    global_op->UpdatePreconditioner();

    // Solve the problem.
    ParameterList lop_list = plist.sublist("solvers").sublist("silent").sublist("pcg parameters");
    AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
        solver(global_op, global_op);
    solver.Init(lop_list);

    CompositeVector& rhs = *global_op->rhs();
    solver.ApplyInverse(rhs, E);

    double heat = op_mag->CalculateOhmicHeating(E);
    double energy = op_mag->CalculateMagneticEnergy(B);
    op_mag->ModifyFields(E, B, dt);

    CHECK(heat > 0.0);
    if (!convergence) CHECK(energy < energy0);
    energy0 = energy;

    cycle++;
    told = tnew;
    tnew += dt;

    // reconstruction 
    B.ScatterMasterToGhosted("face");
    E.ScatterMasterToGhosted("edge");

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

    WhetStone::MFD3D_Electromagnetics mfd(plist, mesh); 
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

    ana.GlobalOp("sum", &avgB, 1);
    ana.GlobalOp("sum", &divB, 1);
    ana.GlobalOp("sum", &errB, 1);

    if (cycle == 1) divB0 = divB;
    CHECK_CLOSE(divB0, divB, 1e-8);

    if (MyPID == 0) {
      std::cout << "time: " << told << "  ||r||=" << solver.residual() 
                << " itr=" << solver.num_itrs() << "  energy= " << energy 
                << "  heat= " << heat
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
      printf("L2(e)=%12.9f  Inf(e)=%11.8f  L2(b)=%12.9f  Inf(b)=%11.8f\n",
          el2_err, einf_err, bl2_err, binf_err);
    }
  }
}

TEST(MAGNETIC_DIFFUSION3D_RELAX) {
  MagneticDiffusion3D<AnalyticElectromagnetics04>(0.1, 0.5, false, 10,8,22, -4.0,-4.0,-10.0, 4.0,4.0,10.0, "structured");
}

TEST(MAGNETIC_DIFFUSION3D_CONVERGENCE) {
  // MagneticDiffusion3D<AnalyticElectromagnetics05>(0.01, 0.1, true, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "structured");
  MagneticDiffusion3D<AnalyticElectromagnetics05>(0.01, 0.1, true, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "test/kershaw08.exo");
}

