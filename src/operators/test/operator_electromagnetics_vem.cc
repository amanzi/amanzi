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
void CurlCurl_VEM(int nx, const std::string& method, int type, int order,
                  double tolerance, const std::string filename = "") {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Curl-curl operator, method=" << method
                            << "  order=" << order << "  type=" << type << std::endl;

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
    // DeformMesh(mesh_tmp, 5, 1.0);
  } else {
    mesh = meshfactory.create(filename, request_faces, request_edges);
  }

  // create resistivity coefficient
  double time = 1.0;
  Analytic ana(1.0, mesh);
  WhetStone::Tensor Kc(3, 2);

  Teuchos::RCP<std::vector<WhetStone::Tensor> > K = Teuchos::rcp(new std::vector<WhetStone::Tensor>());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nedges_wghost = mesh->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    Kc = ana.Tensor(xc, time);
    K->push_back(Kc);
  }

  // create boundary data
  int nde = order + 1;
  int ndf = WhetStone::PolynomialSpaceDimension(2, order) - 1;

  // -- on boundary edges
  Teuchos::RCP<BCs> bc_e = Teuchos::rcp(new BCs(mesh, AmanziMesh::EDGE, WhetStone::DOF_Type::VECTOR));
  auto& bce_model = bc_e->bc_model();
  auto& bce_value = bc_e->bc_value_vector(nde);

  WhetStone::NumericalIntegration numi(mesh);

  for (int e = 0; e < nedges_wghost; ++e) {
    double len = mesh->edge_length(e);
    const auto& tau = mesh->edge_vector(e);
    const auto& xe = mesh->edge_centroid(e);

    ana.set_parameters(tau / len, 0, 0.0);

    if (fabs(xe[0]) < 1e-6 || fabs(xe[0] - 1.0) < 1e-6 ||
        fabs(xe[1]) < 1e-6 || fabs(xe[1] - 1.0) < 1e-6 ||
        fabs(xe[2]) < 1e-6 || fabs(xe[2] - 1.0) < 1e-6) {
      std::vector<double> moments;
      numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 2);

      for (int k = 0; k < moments.size(); ++k) bce_value[e][k] = moments[k];
      bce_model[e] = OPERATOR_BC_DIRICHLET;
    }
  }

  // -- on boundary faces
  Teuchos::RCP<BCs> bc_f;
  if (type == 1 && order == 1) {
    bc_f = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
    auto& bcf_model = bc_f->bc_model();
    auto& bcf_value = bc_f->bc_value_vector(ndf);

    for (int f = 0; f < nfaces_wghost; ++f) {
      double area = mesh->face_area(f);
      const auto& xf = mesh->face_centroid(f);

      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
          fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
          fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) {
        std::vector<double> moments;
        const auto& normal = mesh->face_normal(f);

        WhetStone::SurfaceCoordinateSystem coordsys(xf, normal);

        WhetStone::Polynomial pf(2, order);
        for (auto it = pf.begin(1); it < pf.end(); ++it) {
          int k = it.PolynomialPosition() - 1;
          double factor = std::pow(area, -(double)it.MonomialSetOrder() / 2);
          WhetStone::Polynomial fmono(2, it.multi_index(), factor);
          auto rot = Rot2D(fmono);

          AmanziGeometry::Point tmp = rot[0](0) * (*coordsys.tau())[0] 
                                    + rot[1](0) * (*coordsys.tau())[1];
          ana.set_parameters(tmp, 0, 0.0);
          numi.CalculateFunctionMomentsFace(f, &ana, order - 1, moments, 4);

          bcf_value[f][k] = moments[0];
          bcf_model[f] = OPERATOR_BC_DIRICHLET;
        }
      }
    }
  }

  // create electromagnetics operator
  Teuchos::ParameterList olist = plist->sublist("PK operator").sublist("curlcurl operator");
  olist.sublist("schema").set<std::string>("method", method)
                         .set<int>("method order", order)
                         .set<int>("type", type);
  auto op_curlcurl = Teuchos::rcp(new PDE_Abstract(olist, mesh));
  op_curlcurl->SetBCs(bc_e, bc_e);
  if (type == 1) op_curlcurl->AddBCs(bc_f, bc_f);
  op_curlcurl->Setup(K, false);
  op_curlcurl->UpdateMatrices();

  // add an accumulation operator
  Teuchos::RCP<Operator> global_op = op_curlcurl->global_operator();
  olist = plist->sublist("PK operator").sublist("accumulation operator");
  olist.sublist("schema").set<std::string>("method", method)
                         .set<int>("method order", order)
                         .set<int>("type", type);
  auto op_acc = Teuchos::rcp(new PDE_Abstract(olist, global_op));
  op_acc->SetBCs(bc_e, bc_e);
  if (type == 1) op_acc->AddBCs(bc_f, bc_f);
  op_acc->UpdateMatrices();

  // create source for a manufactured solution.
  const CompositeVectorSpace& cvs = op_curlcurl->global_operator()->DomainMap();
  CompositeVector source(cvs);
  Epetra_MultiVector& src_e = *source.ViewComponent("edge");
  source.PutScalarMasterAndGhosted(0.0);

  std::vector<double> moments;
  for (int c = 0; c < ncells_owned; c++) {
    const auto& Acell = (op_acc->local_op()->matrices)[c]; 
    int nloc = Acell.NumRows();

    WhetStone::DenseVector ana_loc(nloc), src_loc(nloc);

    const auto& edges = mesh->cell_get_edges(c);
    int nedges = edges.size();

    // edge-based DOFs
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      double len = mesh->edge_length(e);
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);

      ana.set_parameters(tau / len, 1, 0.0);

      numi.CalculateFunctionMomentsEdge(e, &ana, order, moments, 3);
      for (int k = 0; k < moments.size(); ++k) ana_loc(n * nde + k) = moments[k];
    }

    // face-based DOFs
    if (type == 1 && order == 1) {
      const auto& faces = mesh->cell_get_faces(c);
      int nfaces = faces.size();

      WhetStone::Polynomial pf(2, order);

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        double area = mesh->face_area(f);
        const auto& xf = mesh->face_centroid(f);
        const auto& normal = mesh->face_normal(f);
        WhetStone::SurfaceCoordinateSystem coordsys(xf, normal);

        for (auto it = pf.begin(1); it < pf.end(); ++it) {
          int k = it.PolynomialPosition() - 1;
          double factor = std::pow(area, -(double)it.MonomialSetOrder() / 2);
          WhetStone::Polynomial fmono(2, it.multi_index(), factor);
          auto rot = Rot2D(fmono);

          AmanziGeometry::Point tmp = rot[0](0) * (*coordsys.tau())[0] 
                                    + rot[1](0) * (*coordsys.tau())[1];
          ana.set_parameters(tmp, 1, 0.0);
          numi.CalculateFunctionMomentsFace(f, &ana, order - 1, moments, 4);
          ana_loc(nedges * nde + ndf * n + k) = moments[0];
        }
      }
    }

    Acell.Multiply(ana_loc, src_loc, false);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      for (int k = 0; k < nde; ++k)
        src_e[k][e] += src_loc(n * nde + k);
    }
    if (type == 1) {
      Epetra_MultiVector& src_f = *source.ViewComponent("face");

      const auto& faces = mesh->cell_get_faces(c);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        for (int k = 0; k < ndf; ++k) 
          src_f[k][f] += src_loc(nedges * nde + n * ndf + k);
      }
    }
  }

  source.GatherGhostedToMaster();

  // set up initial guess for a time-dependent problem
  CompositeVector solution(cvs);
  solution.PutScalar(0.0);

  // BCs, sources, and assemble
  global_op->UpdateRHS(source, true);
  op_curlcurl->ApplyBCs(true, true, true);
  op_acc->ApplyBCs(true, true, false);

  // create preconditioner and test it
  global_op->set_inverse_parameters("Hypre AMG", plist->sublist("preconditioners"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  VerificationCV ver(global_op);
  // ver.CheckMatrixSPD(true, true);
  // ver.CheckPreconditionerSPD(1e-10, true, true);
  // ver.CheckSpectralBounds(0);

  // Create a solver and solve the problem
  CompositeVector& rhs = *global_op->rhs();
  global_op->set_inverse_parameters("Hypre AMG", plist->sublist("preconditioners"), "default", plist->sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ApplyInverse(rhs, solution);

  ver.CheckResidual(solution, 1.0e-10);

  int num_itrs = global_op->num_itrs();
  CHECK(num_itrs < 2000);

  if (MyPID == 0) {
    std::cout << "electric solver: ||r||=" << global_op->residual() 
              << " itr=" << num_itrs
              << " code=" << global_op->returned_code() 
              << " edges=" << rhs.ViewComponent("edge")->GlobalLength() << std::endl;
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
  CurlCurl_VEM<AnalyticElectromagnetics02>(0, "electromagnetics", 2, 0, 1e+2, "test/hexes4.exo");
  CurlCurl_VEM<AnalyticElectromagnetics02>(0, "Nedelec serendipity", 2, 1, 1e+2, "test/hexes4.exo");
  CurlCurl_VEM<AnalyticElectromagnetics02>(0, "Nedelec serendipity", 1, 1, 1e+2, "test/hexes4.exo");
/*
L2(e)= 0.030491648  Inf(e)= 0.033814076  itr= 89  size=10072
L2(e)= 0.015556440  Inf(e)= 0.018429069  itr=162  size=78128
L2(e)= 0.007872589  Inf(e)= 0.008978115  itr=291  size=615520

L2(e)= 0.004197785  Inf(e)= 0.004418447  itr=131  size=10072
L2(e)= 0.001097908  Inf(e)= 0.001145579  itr=243  size=78128
L2(e)= 0.000289011  Inf(e)= 0.000289035  itr=456  size=615520
*/
}

