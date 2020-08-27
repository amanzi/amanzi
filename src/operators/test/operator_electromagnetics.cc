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
#include "MeshFactory.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Electromagnetics.hh"

#include "AnalyticElectromagnetics01.hh"
#include "AnalyticElectromagnetics02.hh"
#include "AnalyticElectromagnetics03.hh"
#include "Verification.hh"

/* *****************************************************************
* TBW 
* **************************************************************** */
template<class Analytic>
void CurlCurl(double c_t, int nx, double tolerance, bool initial_guess,
              const std::string& disc_method = "electromagnetics") {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Curl-curl operator, tol=" << tolerance 
                            << "  method=" << disc_method << std::endl;

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
  if (nx > 0) 
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, request_faces, request_edges);
  else
    mesh = meshfactory.create("test/hex_split_faces5.exo", request_faces, request_edges);

  // create resistivity coefficient
  double time = 1.0;
  Analytic ana(0.0, mesh);
  WhetStone::Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, time);
    K->push_back(Kc);
  }

  // create boundary data
  int nedges_owned = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

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
        bc_value[e] = (ana.electric_exact(xe, time) * tau) / len;
      }
    }
  }

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist.sublist("PK operator").sublist("electromagnetics operator");
  olist.sublist("schema electric").set<std::string>("method", disc_method);
  Teuchos::RCP<PDE_Electromagnetics> op_curlcurl = Teuchos::rcp(new PDE_Electromagnetics(olist, mesh));
  op_curlcurl->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op_curlcurl->global_operator()->DomainMap();

  // create source for a manufactured solution.
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("edge");
  source.PutScalarMasterAndGhosted(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh->cell_get_edges(c, &edges);
    int nedges = edges.size();
    double vol = 3.0 * mesh->cell_volume(c) / nedges;

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      double len = mesh->edge_length(e);
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);
      const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

      src[0][e] += (ana.source_exact(xe, time) * tau) / len * vol;
    }
  }
  source.GatherGhostedToMaster("edge");

  // set up initial guess for a time-dependent problem
  CompositeVector solution(cvs);
  Epetra_MultiVector& sol = *solution.ViewComponent("edge");

  sol.PutScalar(0.0);
  if (initial_guess) {
    for (int e = 0; e < nedges_owned; e++) {
      double len = mesh->edge_length(e);
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);
      const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

      sol[0][e] = (ana.electric_exact(xe, time) * tau) / len;
    }
  } 

  // set up the diffusion operator
  op_curlcurl->SetTensorCoefficient(K);
  op_curlcurl->UpdateMatrices();

  // Add an accumulation term.
  CompositeVector phi(cvs);
  phi.PutScalar(c_t);

  Teuchos::RCP<Operator> global_op = op_curlcurl->global_operator();
  Teuchos::RCP<PDE_Accumulation> op_acc = Teuchos::rcp(new PDE_Accumulation(AmanziMesh::EDGE, global_op));

  double dT = 1.0;
  op_acc->AddAccumulationDelta(solution, phi, phi, dT, "edge");

  // BCs, sources, and assemble
  op_curlcurl->ApplyBCs(true, true, true);
  global_op->UpdateRHS(source, false);

  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);
  ver.CheckPreconditionerSPD(1e-12, true, true);

  // re init with a solver...
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "default", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();
  
  CompositeVector& rhs = *global_op->rhs();
  global_op->ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-10);

  int num_itrs = global_op->num_itrs();
  CHECK(num_itrs < 100);

  if (MyPID == 0) {
    std::cout << "electric solver (pcg): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() << std::endl;
  }

  // compute electric error
  Epetra_MultiVector& E = *solution.ViewComponent("edge", true);
  double enorm, el2_err, einf_err;
  ana.ComputeEdgeError(E, time, enorm, el2_err, einf_err);

  if (MyPID == 0) {
    el2_err /= enorm;
    printf("L2(e)=%12.9f  Inf(e)=%9.6f  itr=%3d  size=%d\n", el2_err, einf_err,
            global_op->num_itrs(), rhs.GlobalLength());

    CHECK(el2_err < tolerance);
  }
}


TEST(CURL_CURL_LINEAR) {
  CurlCurl<AnalyticElectromagnetics01>(1.0e-5, 0, 1e-4, false);
  // CurlCurl<AnalyticElectromagnetics01>(1.0e-5, 0, 1e-4, false, "mfd: generalized");
}

/*
TEST(CURL_CURL_NONLINEAR) {
  CurlCurl<AnalyticElectromagnetics02>(1.0e-1, 0, 2e-1, false);
}

TEST(CURL_CURL_TIME_DEPENDENT) {
  CurlCurl<AnalyticElectromagnetics03>(1.0, 0, 2e-3, true);
  // CurlCurl<AnalyticElectromagnetics03>(1.0, 0, 2e-3, true, "mfd: generalized");
}
*/

