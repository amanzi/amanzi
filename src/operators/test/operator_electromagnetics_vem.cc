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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "WhetStoneFunction.hh"

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
void CurlCurl_VEM(int nx, const std::string& method, int order, double tolerance, 
                  const std::string filename = "") {
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
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a MSTK mesh framework
  ParameterList region_list = plist->sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory meshfactory(comm, gm, mesh_list);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  bool request_faces(true), request_edges(true);
  RCP<const Mesh> mesh;
  if (nx > 0) {
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, request_faces, request_edges);
    auto mesh_tmp = Teuchos::rcp_const_cast<Mesh>(mesh);
    DeformMesh(mesh_tmp, 5, 1.0);
  } else {
    mesh = meshfactory.create(filename, request_faces, request_edges);
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
  int nde = order + 1;
  int nedges_wghost = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, WhetStone::DOF_Type::VECTOR));
  auto& bc_model = bc->bc_model();
  auto& bc_value = bc->bc_value_vector(nde);

  std::vector<int> edirs;
  AmanziMesh::Entity_ID_List cells, edges;
  WhetStone::NumericalIntegration<AmanziMesh::Mesh> numi(mesh);

  for (int e = 0; e < nedges_wghost; ++e) {
    double len = mesh->edge_length(e);
    const auto& tau = mesh->edge_vector(e);
    const auto& xe = mesh->edge_centroid(e);

    std::vector<AmanziGeometry::Point> coordsys(1, tau);
    ana.set_parameters(tau / len, 0, 0.0);

    if (fabs(xe[0]) < 1e-6 || fabs(xe[0] - 1.0) < 1e-6 ||
        fabs(xe[1]) < 1e-6 || fabs(xe[1] - 1.0) < 1e-6 ||
        fabs(xe[2]) < 1e-6 || fabs(xe[2] - 1.0) < 1e-6) {
      std::vector<double> moments;
      numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 2);

      for (int k = 0; k < moments.size(); ++k) bc_value[e][k] = moments[k];
      bc_model[e] = OPERATOR_BC_DIRICHLET;
    }
  }

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist->sublist("PK operator").sublist("curlcurl operator");
  olist.sublist("schema").set<std::string>("method", method);
  olist.sublist("schema").set<int>("method order", order);
  auto op_curlcurl = Teuchos::rcp(new PDE_Abstract(olist, mesh));
  op_curlcurl->SetBCs(bc, bc);
  op_curlcurl->Setup(K, false);
  op_curlcurl->UpdateMatrices();

  // add an accumulation operator
  Teuchos::RCP<Operator> global_op = op_curlcurl->global_operator();
  olist = plist->sublist("PK operator").sublist("accumulation operator");
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
    int nloc = nedges * nde;
    WhetStone::DenseVector ana_loc(nloc, 0.0), src_loc(nloc);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      double len = mesh->edge_length(e);
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);

      std::vector<AmanziGeometry::Point> coordsys(1, tau);
      ana.set_parameters(tau / len, 1, 0.0);

      std::vector<double> moments;
      numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 3);
      for (int k = 0; k < moments.size(); ++k) ana_loc(n * nde + k) = moments[k];
    }

    const auto& Acell = (op_acc->local_op()->matrices)[c]; 
    Acell.Multiply(ana_loc, src_loc, false);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      for (int k = 0; k < nde; ++k) 
        src[k][e] += src_loc(n * nde + k);
    }
  }

  source.GatherGhostedToMaster("edge");

  // set up initial guess for a time-dependent problem
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // BCs, sources, and assemble
  global_op->UpdateRHS(source, true);
  op_curlcurl->ApplyBCs(true, true, true);
  op_acc->ApplyBCs(true, true, false);

  // create preconditoner
  global_op->set_inverse_parameters("Hypre AMG", plist->sublist("preconditioners"), 
                                    "silent", plist->sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);
  if (comm->NumProc() == 1)
    ver.CheckPreconditionerSPD(1e-10, true, true);

  // Solve the problem.
  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-10);

  int num_itrs = global_op->num_itrs();
  CHECK(num_itrs < 2000);

  if (MyPID == 0) {
    std::cout << "electric solver: ||r||=" << global_op->residual() 
              << " itr=" << num_itrs
              << " code=" << global_op->returned_code() << std::endl;
  }

  // compute electric error
  Epetra_MultiVector& E = *solution.ViewComponent("edge", true);
  double enorm, el2_err, einf_err;
  ana.ComputeEdgeError(E, time, enorm, el2_err, einf_err);

  if (MyPID == 0) {
    el2_err /= enorm;
    printf("L2(e)=%12.9f  Inf(e)=%12.9f  itr=%3d  size=%d\n",
            el2_err, einf_err, global_op->num_itrs(), rhs.GlobalLength() * nde);

    CHECK(el2_err < tolerance);
  }
}


TEST(CURL_CURL_HIGH_ORDER) {
  // CurlCurl_VEM<AnalyticElectromagnetics02>(0, "electromagnetics", 0, 1e+3, "test/hexes8.exo");
  CurlCurl_VEM<AnalyticElectromagnetics02>(0, "Nedelec serendipity type2", 1, 1e+2, "test/hexes4.exo");
  // CurlCurl_VEM<AnalyticElectromagnetics02>(0, "Nedelec serendipity type2", 1, 1e+2, "test/hexes16.exo");
/*
L2(e)= 0.030491648  Inf(e)= 0.033814076  itr= 89  size=10072
L2(e)= 0.015556440  Inf(e)= 0.018429069  itr=162  size=78128
L2(e)= 0.007872589  Inf(e)= 0.008978115  itr=291  size=615520

L2(e)= 0.006310242  Inf(e)= 0.006942832  itr=249  size=20144
L2(e)= 0.001723041  Inf(e)= 0.001851226  itr=401  size=156256

L2(e)= 0.004197785  Inf(e)= 0.004418447  itr=131  size=10072
L2(e)= 0.001097908  Inf(e)= 0.001145579  itr=243  size=78128
L2(e)= 0.000289011  Inf(e)= 0.000289035  itr=456  size=615520
  CurlCurl_VEM<AnalyticElectromagnetics01>(0, "electromagnetics", 0, 1e-8, "test/hex_split_faces5.exo");
  CurlCurl_VEM<AnalyticElectromagnetics02>(8, "electromagnetics", 0, 1e-3);
  CurlCurl_VEM<AnalyticElectromagnetics02>(7, "Nedelec serendipity type2", 0, 2e-3);
*/
}

