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
#include "MeshFactory.hh"
#include "MFD3D_Electromagnetics.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStoneFunction.hh"
#include "WhetStoneMeshUtils.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_MagneticDiffusion.hh"

#include "AnalyticElectromagnetics04.hh"
#include "AnalyticElectromagnetics05.hh"
#include "MeshDeformation.hh"

/* *****************************************************************
* Testing operators for Maxwell-type problems: 3D
* **************************************************************** */
template<class Analytic>
void MagneticDiffusionVEM(
    double dt, double tend, int order,
    int nx, int ny, int nz,
    double Xa, double Ya, double Za, double Xb, double Yb, double Zb,
    const std::string& name)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::WhetStone;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Magnetic diffusion VEM, dt=" 
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
  if (name == "structured") {
    mesh = meshfactory.create(Xa, Ya, Za, Xb, Yb, Zb, nx, ny, nz, request_faces, request_edges);
    // DeformMesh(mesh, 5, 1.0);
  } else {
    mesh = meshfactory.create(name, request_faces, request_edges);
    // mesh = meshfactory.create("test/hex_split_faces5.exo", request_faces, request_edges);
  }

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nedges_owned = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Analytic ana(mesh);

  Polynomial pe(1, order), pf(2, order);
  int nde = PolynomialSpaceDimension(1, order);
  int ndf = PolynomialSpaceDimension(2, order);

  // create resistivity coefficient
  double told(0.0), tnew(dt);
  Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<Tensor> > K = Teuchos::rcp(new std::vector<Tensor>());

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, tnew);
    K->push_back(Kc);
  }

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("magnetic diffusion operators vem");
  olist.sublist("schema electric").set<int>("method order", order);
  olist.sublist("schema magnetic").set<int>("method order", order);

  Teuchos::RCP<PDE_MagneticDiffusion> op_mag = Teuchos::rcp(new PDE_MagneticDiffusion(olist, mesh));
  Teuchos::RCP<Operator> global_op = op_mag->global_operator();

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, WhetStone::DOF_Type::VECTOR));
  op_mag->SetBCs(bc, bc);

  // create/extract solution maps
  const CompositeVectorSpace& cvs_e = global_op->DomainMap();
  Teuchos::RCP<CompositeVectorSpace> cvs_b = Teuchos::rcp(new CompositeVectorSpace());
  cvs_b->SetMesh(mesh)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, ndf);

  CompositeVector E(cvs_e);
  CompositeVector B(*cvs_b);

  // set up initial guess for a time-dependent problem 
  WhetStone::NumericalIntegration<AmanziMesh::Mesh> numi(mesh);

  Epetra_MultiVector& Ee = *E.ViewComponent("edge");
  Epetra_MultiVector& Bf = *B.ViewComponent("face");

  Ee.PutScalar(0.0);
  Bf.PutScalar(0.0);

  std::vector<const WhetStone::WhetStoneFunction*> funcs(2);
  funcs[0] = &ana;

  for (int e = 0; e < nedges_owned; ++e) {
    double len = mesh->edge_length(e);
    const AmanziGeometry::Point& tau = mesh->edge_vector(e);
    const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

    std::vector<AmanziGeometry::Point> coordsys(1, tau);
    ana.set_parameters(tau / len, 0, told);

    for (auto it = pe.begin(); it < pe.end(); ++it) {
      const int* index = it.multi_index();
      Polynomial emono(1, index, 1.0);
      emono.InverseChangeCoordinates(xe, coordsys);  
      funcs[1] = &emono;

      int k = it.PolynomialPosition();
      Ee[k][e] = numi.IntegrateFunctionsEdge(e, funcs, 2) / len;
    }
  }

  for (int f = 0; f < nfaces_owned; ++f) {
    double area = mesh->face_area(f);
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
 
    SurfaceCoordinateSystem coordsys(xf, normal);
    ana.set_parameters(normal / area, 2, told);

    for (auto it = pf.begin(); it < pf.end(); ++it) {
      int m = it.MonomialSetOrder();
      double scale = std::pow(area, -(double)m / 2);

      const int* index = it.multi_index();
      Polynomial fmono(2, index, scale);
      fmono.InverseChangeCoordinates(xf, *coordsys.tau());  

      funcs[1] = &fmono;

      int k = it.PolynomialPosition();
      Bf[k][f] = numi.IntegrateFunctionsTriangulatedFace(f, funcs, 1) / area;
    }
  }

  int cycle(0);
  double energy0(1e+99), divB0(0.0);
  while (told + dt/2 < tend) {
    // set up the diffusion operator
    global_op->Init();
    op_mag->SetTensorCoefficient(K);
    op_mag->UpdateMatrices();

    // update BCs
    auto& bc_model = bc->bc_model();
    auto& bc_value = bc->bc_value_vector(nde);

    std::vector<int> edirs;
    AmanziMesh::Entity_ID_List cells, edges;

    for (int f = 0; f < nfaces_wghost; ++f) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);

      if (fabs(xf[0] - Xa) < 1e-6 || fabs(xf[0] - Xb) < 1e-6 ||
          fabs(xf[1] - Ya) < 1e-6 || fabs(xf[1] - Yb) < 1e-6 ||
          fabs(xf[2] - Za) < 1e-6 || fabs(xf[2] - Zb) < 1e-6) {
        mesh->face_get_edges_and_dirs(f, &edges, &edirs);
        int nedges = edges.size();
        for (int i = 0; i < nedges; ++i) {
          int e = edges[i];
          double len = mesh->edge_length(e);
          const AmanziGeometry::Point& tau = mesh->edge_vector(e);
          const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

          std::vector<AmanziGeometry::Point> coordsys(1, tau);
          ana.set_parameters(tau / len, 0, tnew - dt/2);

          for (auto it = pe.begin(); it < pe.end(); ++it) {
            const int* index = it.multi_index();
            Polynomial emono(1, index, 1.0);
            emono.InverseChangeCoordinates(xe, coordsys);  
            funcs[1] = &emono;

            int k = it.PolynomialPosition();
            bc_value[e][k] = numi.IntegrateFunctionsEdge(e, funcs, 2) / len;
          }
          bc_model[e] = OPERATOR_BC_DIRICHLET;
        }
      }
    }

    // BCs, sources, and assemble
    op_mag->ModifyMatrices(E, B, dt);
    op_mag->ApplyBCs(true, true, true);

    // Solve the problem.
    global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "silent", plist.sublist("solvers"));
    global_op->InitializeInverse();
    global_op->ComputeInverse();

    CompositeVector& rhs = *global_op->rhs();
    global_op->ApplyInverse(rhs, E);

    double heat = op_mag->CalculateOhmicHeating(E);
    double energy = op_mag->CalculateMagneticEnergy(B);
    op_mag->ModifyFields(E, B, dt);

    CHECK(heat > 0.0);
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
      std::cout << "time: " << told << "  ||r||=" << global_op->residual() 
                << " itr=" << global_op->num_itrs() << "  energy= " << energy 
                << "  heat= " << heat
                << "  avgB=" << avgB / ncells_owned 
                << "  divB=" << std::pow(divB, 0.5) 
                << "  errB=" << std::pow(errB, 0.5) << std::endl;
    }
  }

  // compute electric and magnetic errors
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

TEST(MAGNETIC_DIFFUSION3D_CONVERGENCE) {
  MagneticDiffusionVEM<AnalyticElectromagnetics05>(0.01, 0.1, 0, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "structured");
  // MagneticDiffusionVEM<AnalyticElectromagnetics05>(0.01, 0.1, 0, 8,8,8, 0.0,0.0,0.0, 1.0,1.0,1.0, "test/kershaw08.exo");
}

