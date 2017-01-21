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
#include "LinearOperatorFactory.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mfd3d_diffusion.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Accumulation.hh"
#include "ElectromagneticsMHD.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

#include "AnalyticElectromagnetics04.hh"

/* *****************************************************************
* TBW 
* **************************************************************** */
template<class Analytic>
void ResistiveMHD(double dt, double tend, bool initial_guess) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Resistive MHD" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_electromagnetics.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a MSTK mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  double Xb(4.0), Yb(4.0), Zb(10.0);
  bool request_faces(true), request_edges(true);
  RCP<const Mesh> mesh = meshfactory(-Xb, -Yb, -Zb, Xb, Yb, Zb, 10, 8, 22, gm, request_faces, request_edges);
  // RCP<const Mesh> mesh = meshfactory("test/hex_split_faces5.exo", gm, request_faces, request_edges);
  // RCP<const Mesh> mesh = meshfactory("test/isohelix.exo", gm, request_faces, request_edges);

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

  Teuchos::RCP<BCs> bc1 = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE));
  std::vector<int>& bc_model = bc1->bc_model();
  std::vector<double>& bc_value = bc1->bc_value();

  Teuchos::RCP<BCs> bc2 = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE));
  std::vector<int>& bc_model2 = bc2->bc_model();
  std::vector<Point>& bc_value2 = bc2->bc_value_point();

  std::vector<int> edirs;
  AmanziMesh::Entity_ID_List cells, edges;

  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0] + Xb) < 1e-6) {
      bc_model2[f] = OPERATOR_BC_NEUMANN;
      bc_value2[f] = Point(0.0, 0.0, 1.0);
    } 
    else if (fabs(xf[0] - Xb) < 1e-6 ||
    // if (fabs(xf[0] + Xb) < 1e-6 || fabs(xf[0] - Xb) < 1e-6 ||
            fabs(xf[1] + Yb) < 1e-6 || fabs(xf[1] - Yb) < 1e-6 ||
            fabs(xf[2] + Zb) < 1e-6 || fabs(xf[2] - Zb) < 1e-6) {
      mesh->face_get_edges_and_dirs(f, &edges, &edirs);
      int nedges = edges.size();
      for (int i = 0; i < nedges; ++i) {
        int e = edges[i];
        double len = mesh->edge_length(e);
        const AmanziGeometry::Point& tau = mesh->edge_vector(e);
        const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

        bc_model[e] = OPERATOR_BC_DIRICHLET;
        bc_value[e] = (ana.electric_exact(xe, tnew) * tau) / len;
      }
    }
  }

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("electromagnetics operator");
  Teuchos::RCP<ElectromagneticsMHD> op_mhd = Teuchos::rcp(new ElectromagneticsMHD(olist, mesh));
  op_mhd->SetBCs(bc1, bc1);
  op_mhd->AddBCs(bc2, bc2);

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

  if (initial_guess) {
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
  } 
  // CompositeVector B0(B);
  // Epetra_MultiVector& B0f = *B0.ViewComponent("face");

  int cycle(0);
  double energy0(1e+99), divB0(0.0);
  while (told < tend) {
    // set up the diffusion operator
    global_op->Init();
    op_mhd->SetTensorCoefficient(K);
    op_mhd->UpdateMatrices();

    // Add an accumulation term using dt=1 since time step is taken into
    // account in the system modification routine.
    CompositeVector phi(cvs_e);
    phi.PutScalar(1.0);

    Teuchos::RCP<Accumulation> op_acc =
        Teuchos::rcp(new Accumulation(AmanziMesh::EDGE, global_op));

    op_acc->AddAccumulationTerm(phi, 1.0, "edge");

    // BCs, sources, and assemble
    op_mhd->ModifyMatrices(E, B, dt);
    op_mhd->ApplyBCs(true, true);
    op_acc->ApplyBCs(bc1);
    global_op->SymbolicAssembleMatrix();
    global_op->AssembleMatrix();

    ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
    global_op->InitPreconditioner("Hypre AMG", slist);

    // Solve the problem.
    ParameterList lop_list = plist.get<Teuchos::ParameterList>("solvers");
    AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
        solver = factory.Create("silent", lop_list, global_op);

    CompositeVector& rhs = *global_op->rhs();
    int ierr = solver->ApplyInverse(rhs, E);

    double energy = op_mhd->CalculateMagneticEnergy(B);
    op_mhd->ModifyFields(E, B, dt);

    CHECK(energy < energy0);
    energy0 = energy;

    cycle++;
    told = tnew;
    tnew += dt;

    // reconstruction
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 3);

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
        for (int k = 0; k < 3; ++k) {
          sol[k][c] += Bf[0][f] * dirs[n] * area * (xf[k] - xc[k]) / vol;
        }        
        tmp += Bf[0][f] * dirs[n] * area / vol;
      }
      avgB += std::fabs(sol[0][c]);
      divB += tmp * tmp * vol; 
      errB += vol * (std::pow(sol[0][c], 2.0) + std::pow(sol[1][c], 2.0) 
                   + std::pow(sol[2][c] - 1.0, 2.0));
    }

    if (cycle == 1) divB0 = divB;
    CHECK_CLOSE(divB0, divB, 1e-8);

    if (MyPID == 0) {
      std::cout << "time: " << told << "  ||r||=" << solver->residual() 
                << " itr=" << solver->num_itrs() << "  energy= " << energy 
                << "  avgB=" << avgB / ncells_owned 
                << "  divB=" << std::pow(divB, 0.5) 
                << "  ||B-1||=" << std::pow(errB, 0.5) << std::endl;
    }

    // visualization
    if (MyPID == 0 && (cycle % 5 == 0)) {
      GMV::open_data_file(*mesh, "operators.gmv");
      GMV::start_data();
      GMV::write_cell_data(sol, 0, "Bx");
      GMV::write_cell_data(sol, 1, "By");
      GMV::write_cell_data(sol, 2, "Bz");
      GMV::close_data_file();
    }
  }
}


TEST(RESISTIVE_MHD_LINEAR) {
  ResistiveMHD<AnalyticElectromagnetics04>(0.1, 0.5, true);
}

