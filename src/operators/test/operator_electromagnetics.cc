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
#include "Electromagnetics.hh"
#include "OperatorDefs.hh"

#include "AnalyticElectromagnetics01.hh"
#include "AnalyticElectromagnetics02.hh"
#include "AnalyticElectromagnetics03.hh"
#include "Verification.hh"

/* *****************************************************************
* TBW 
* **************************************************************** */
template<class Analytic>
void CurlCurl(double c_t, double tolerance, bool initial_guess) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Curl-curl operator, tol=" << tolerance << std::endl;

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

  bool request_faces(true), request_edges(true);
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, gm, request_faces, request_edges);
  RCP<const Mesh> mesh = meshfactory("test/hex_split_faces5.exo", gm, request_faces, request_edges);

  // create resistivity coefficient
  double time = 1.0;
  Analytic ana(mesh);
  WhetStone::Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, time);
    K->push_back(Kc);
  }

  // create boundary data
  int nedges_owned = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::OWNED);
  int nedges_wghost = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::USED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE));
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
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("electromagnetics operator");
  Teuchos::RCP<Electromagnetics> op_curlcurl = Teuchos::rcp(new Electromagnetics(olist, mesh));
  op_curlcurl->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op_curlcurl->global_operator()->DomainMap();

  // create source for a manufactured solution.
  CompositeVector source(cvs);
  Epetra_MultiVector& src = *source.ViewComponent("edge");

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
  Teuchos::RCP<Accumulation> op_acc = Teuchos::rcp(new Accumulation(AmanziMesh::EDGE, global_op));

  double dT = 1.0;
  op_acc->AddAccumulationDelta(solution, phi, phi, dT, "edge");

  // BCs, sources, and assemble
  op_curlcurl->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();
  global_op->UpdateRHS(source, false);

  ParameterList slist = plist.get<Teuchos::ParameterList>("preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  // Test SPD properties of the matrix and preconditioner.
  VerificationCV ver(global_op);
  ver.CheckMatrixSPD(true, true);
  ver.CheckPreconditionerSPD(true, true);

  // Solve the problem.
  ParameterList lop_list = plist.get<Teuchos::ParameterList>("solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create("default", lop_list, global_op);

  CompositeVector& rhs = *global_op->rhs();
  int ierr = solver->ApplyInverse(rhs, solution);

  int num_itrs = solver->num_itrs();
  CHECK(num_itrs < 100);

  if (MyPID == 0) {
    std::cout << "electric solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;
  }

  // compute electric error
  Epetra_MultiVector& E = *solution.ViewComponent("edge", false);
  double enorm, el2_err, einf_err;
  ana.ComputeEdgeError(E, time, enorm, el2_err, einf_err);

  if (MyPID == 0) {
    el2_err /= enorm;
    printf("L2(e)=%10.7f  Inf(e)=%9.6f  itr=%3d  size=%d\n", el2_err, einf_err,
            solver->num_itrs(), rhs.GlobalLength());

    CHECK(el2_err < tolerance);
  }
}


TEST(CURL_CURL_LINEAR) {
  CurlCurl<AnalyticElectromagnetics01>(1.0e-3, 1e-3, false);
}

TEST(CURL_CURL_NONLINEAR) {
  CurlCurl<AnalyticElectromagnetics02>(1.0e-3, 2e-1, false);
}

TEST(CURL_CURL_TIME_DEPENDENT) {
  CurlCurl<AnalyticElectromagnetics03>(1.0, 2e-3, true);
}
