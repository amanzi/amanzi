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
#include "Tensor.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"

#include "AnalyticElectromagnetics01.hh"
#include "AnalyticElectromagnetics02.hh"
#include "MeshDeformation.hh"
#include "Verification.hh"

/* *****************************************************************
* TBW 
* **************************************************************** */
template<class Analytic>
void CurlCurl_VEM(int nx, const std::string method, int order, double tolerance) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Curl-curl operator, method=" << method
                            << "  order=" << order << std::endl;

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
  RCP<const Mesh> mesh;
  if (nx > 0) {
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, request_faces, request_edges);
    auto mesh_tmp = Teuchos::rcp_const_cast<Mesh>(mesh);
    DeformMesh(mesh_tmp, 5, 1.0);
  } else {
    mesh = meshfactory.create("test/hex_split_faces5.exo", request_faces, request_edges);
  }

  // create resistivity coefficient
  double time = 1.0;
  Analytic ana(1.0, mesh);
  WhetStone::Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, time);
    K->push_back(Kc);
  }

  // create boundary data
  // int nedges_owned = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, WhetStone::DOF_Type::VECTOR));
  auto& bc_model = bc->bc_model();
  auto& bc_value = bc->bc_value_vector(order + 1);

  std::vector<int> edirs;
  AmanziMesh::Entity_ID_List cells, edges;

  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);

    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {

      mesh->face_get_edges_and_dirs(f, &edges, &edirs);
      int nedges = edges.size();
      for (int i = 0; i < nedges; ++i) {
        int e = edges[i];
        double len = mesh->edge_length(e);
        const AmanziGeometry::Point& tau = mesh->edge_vector(e);
        const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

        bc_model[e] = OPERATOR_BC_DIRICHLET;
        bc_value[e][0] = (ana.electric_exact(xe, time) * tau) / len;
        // if (order > 0) bc_value[e][1] = 0.0;
      }
    }
  }

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("curlcurl operator");
  olist.sublist("schema").set<std::string>("method", method);
  olist.sublist("schema").set<int>("method order", order);
  auto op_curlcurl = Teuchos::rcp(new PDE_Abstract(olist, mesh));
  op_curlcurl->SetBCs(bc, bc);
  op_curlcurl->Setup(K, false);
  op_curlcurl->UpdateMatrices();

  // add an accumulation aoperator
  Teuchos::RCP<Operator> global_op = op_curlcurl->global_operator();
  olist = plist.sublist("PK operator").sublist("accumulation operator");
  olist.sublist("schema").set<std::string>("method", method);
  olist.sublist("schema").set<int>("method order", order);
  auto op_acc = Teuchos::rcp(new PDE_Abstract(olist, global_op));
  op_acc->SetBCs(bc, bc);
  op_acc->UpdateMatrices();

  // create source for a manufactured solution.
  const CompositeVectorSpace& cvs = op_curlcurl->global_operator()->DomainMap();
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("edge");
  source.PutScalarMasterAndGhosted(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh->cell_get_edges(c, &edges);
    int nedges = edges.size();
    WhetStone::DenseVector ana_loc(nedges), src_loc(nedges);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      double len = mesh->edge_length(e);
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);
      const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

      ana_loc(n) = (ana.source_exact(xe, time) * tau) / len;
    }

    const auto& Acell = (op_acc->local_op()->matrices)[c]; 
    Acell.Multiply(ana_loc, src_loc, false);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      src[0][e] += src_loc(n);
      if (order > 0) src[1][e] += 0.0;
    }
  }

  source.GatherGhostedToMaster("edge");

  // set up initial guess for a time-dependent problem
  CompositeVector solution(cvs);
  Epetra_MultiVector& sol = *solution.ViewComponent("edge");
  sol.PutScalar(0.0);

  // BCs, sources, and assemble
  global_op->UpdateRHS(source, true);

  op_curlcurl->ApplyBCs(true, true, true);
  op_acc->ApplyBCs(true, true, false);

  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  ParameterList slist = plist.sublist("preconditioners").sublist("Hypre AMG");
  global_op->InitializePreconditioner(slist);
  global_op->UpdatePreconditioner();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);
  ver.CheckPreconditionerSPD(1e-10, true, true);

  // Solve the problem.
  ParameterList lop_list = plist.sublist("solvers").sublist("default").sublist("pcg parameters");
  AmanziSolvers::LinearOperatorPCG<Operator, CompositeVector, CompositeVectorSpace>
      solver(global_op, global_op);
  solver.Init(lop_list);

  CompositeVector& rhs = *global_op->rhs();
  solver.ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-10);

  int num_itrs = solver.num_itrs();
  CHECK(num_itrs < 100);

  if (MyPID == 0) {
    std::cout << "electric solver (pcg): ||r||=" << solver.residual() 
              << " itr=" << solver.num_itrs()
              << " code=" << solver.returned_code() << std::endl;
  }

  // compute electric error
  Epetra_MultiVector& E = *solution.ViewComponent("edge", true);
  double enorm, el2_err, einf_err;
  ana.ComputeEdgeError(E, time, enorm, el2_err, einf_err);

  if (MyPID == 0) {
    el2_err /= enorm;
    printf("L2(e)=%12.9f  Inf(e)=%12.9f  itr=%3d  size=%d\n",
            el2_err, einf_err, solver.num_itrs(), rhs.GlobalLength());

    CHECK(el2_err < tolerance);
  }
}


TEST(CURL_CURL_HIGH_ORDER) {
  CurlCurl_VEM<AnalyticElectromagnetics01>(0, "electromagnetics", 0, 1e-8);
  CurlCurl_VEM<AnalyticElectromagnetics02>(8, "electromagnetics", 0, 1e-3);
  // CurlCurl_VEM<AnalyticElectromagnetics02>(8, "Nedelec serendipity type2", 0, 2e-3);
  CurlCurl_VEM<AnalyticElectromagnetics02>(8, "Nedelec serendipity type2", 1, 2e-3);
}

