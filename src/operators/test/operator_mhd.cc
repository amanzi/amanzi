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
#include "AnalyticMHD_01.hh"
#include "AnalyticMHD_02.hh"
#include "AnalyticMHD_03.hh"

#include "OperatorAccumulation.hh"
#include "OperatorCurlCurl.hh"
#include "OperatorDefs.hh"
#include "Verification.hh"

/* *****************************************************************
* TBW 
* **************************************************************** */
template<class Analytic>
void ResistiveMHD(double c_t, double tolerance) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Resistive MHD" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_mhd.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a MSTK mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  bool request_faces(true), request_edges(true);
  // RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 8, 8, gm, request_faces, request_edges);
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
  std::vector<int> bc_model(nedges_wghost, OPERATOR_BC_NONE);
  std::vector<double> bc_value(nedges_wghost);
  std::vector<double> bc_mixed;

  for (int e = 0; e < nedges_wghost; e++) {
    double len = mesh->edge_length(e);
    const AmanziGeometry::Point& tau = mesh->edge_vector(e);
    const AmanziGeometry::Point& xe = mesh->edge_centroid(e);

    if (fabs(xe[0]) < 1e-6 || fabs(xe[0] - 1.0) < 1e-6 ||
        fabs(xe[1]) < 1e-6 || fabs(xe[1] - 1.0) < 1e-6 ||
        fabs(xe[2]) < 1e-6 || fabs(xe[2] - 1.0) < 1e-6) {
      bc_model[e] = OPERATOR_BC_DIRICHLET;
      bc_value[e] = (ana.electric_exact(xe, time) * tau) / len;
    }
  }
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(OPERATOR_BC_TYPE_EDGE, bc_model, bc_value, bc_mixed));

  // create curl-curl operator
  Teuchos::ParameterList olist = plist.get<Teuchos::ParameterList>("PK operator")
                                      .get<Teuchos::ParameterList>("curl-curl operator");
  Teuchos::RCP<OperatorCurlCurl> op_curlcurl = Teuchos::rcp(new OperatorCurlCurl(olist, mesh));
  op_curlcurl->SetBCs(bc, bc);
  const CompositeVectorSpace& cvs = op_curlcurl->global_operator()->DomainMap();

  // create source for a manufactured solution.
  AmanziMesh::Entity_ID_List edges;
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
  if (c_t != 0.0) {
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
  Teuchos::RCP<OperatorAccumulation> op_acc =
      Teuchos::rcp(new OperatorAccumulation(AmanziMesh::EDGE, global_op));

  double dT = 1.0;
  op_acc->AddAccumulationTerm(solution, phi, dT, "edge");

  // BCs, sources, and assemble
  op_curlcurl->ApplyBCs(true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();
  global_op->UpdateRHS(source, false);

  // Create a preconditioner.
  ParameterList slist = plist.get<Teuchos::ParameterList>("Preconditioners");
  global_op->InitPreconditioner("Hypre AMG", slist);

  // Test SPD properties of the matrix and preconditioner.
  Verification ver(global_op);
  ver.CheckMatrixSPD(false, true);
  ver.CheckPreconditionerSPD(false, true);

  // Solve the problem.
  ParameterList lop_list = plist.get<Teuchos::ParameterList>("Solvers");
  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create("AztecOO CG", lop_list, global_op);

  CompositeVector& rhs = *global_op->rhs();
  int ierr = solver->ApplyInverse(rhs, solution);

  int num_itrs = solver->num_itrs();
  CHECK(num_itrs < 100);

  if (MyPID == 0) {
    std::cout << "pressure solver (" << solver->name() 
              << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
              << " code=" << solver->returned_code() << std::endl;
  }

  // compute electric error
  Epetra_MultiVector& E = *solution.ViewComponent("edge", false);
  double enorm, el2_err, einf_err;
  ana.ComputeEdgeError(E, time, enorm, el2_err, einf_err);

  if (MyPID == 0) {
    el2_err /= enorm; 
    el2_err /= enorm;
    printf("L2(e)=%9.6f  Inf(e)=%9.6f  itr=%3d  size=%d\n", el2_err, einf_err,
            solver->num_itrs(), rhs.GlobalLength());

    CHECK(el2_err < tolerance);
  }
}

TEST(RESISTIVE_MHD_LINEAR) {
  ResistiveMHD<AnalyticMHD_01>(0.0, 1e-2);
}

TEST(RESISTIVE_MHD_NONLINEAR) {
  ResistiveMHD<AnalyticMHD_02>(0.0, 2e-2);
}

TEST(RESISTIVE_MHD_TIME_DEPENDENT) {
  ResistiveMHD<AnalyticMHD_03>(1.0, 2e-2);
}
