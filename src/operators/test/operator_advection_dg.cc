/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  DG methods for linear advection equations.
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
#include "Explicit_TI_RK.hh"
#include "GMVMesh.hh"
#include "CompositeVector.hh"
#include "DG_Modal.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "AnalyticDG02.hh"
#include "AnalyticDG03.hh"

#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"


/* *****************************************************************
* This tests exactness of the advection scheme for steady-state
* equations c p + div(v p) = f  (conservative formulation)
* and       c p + v . grad(p) = f.
* Two ways to impose Dirichlet BCs are used: primal and dual.
* **************************************************************** */
template<class AnalyticDG>
void AdvectionSteady(int dim, std::string filename, int nx,
                     std::string weak_form, bool conservative_form,
                     std::string dg_basis = "regularized")
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  std::string problem = (conservative_form) ? ", conservative formulation" : "";
  if (MyPID == 0) std::cout << "\nTest: " << dim 
                            << "D steady advection, dG method" << problem
                            << ", weak formulation=" << weak_form
                            << ", basis=" << dg_basis << std::endl;

  // read parameter list
  std::string xmlFileName;
  xmlFileName = "test/operator_advection_dg.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create a mesh framework
  Teuchos::RCP<GeometricModel> gm;
  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  double weak_sign = 1.0;
  std::string pk_name;
  RCP<const Mesh> mesh;

  if (dim == 2) {
    // mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, nx);
    mesh = meshfactory.create(filename, true, true);
    pk_name = "PK operator 2D";

    if (weak_form == "primal") {
      weak_sign = -1.0;
      pk_name = "PK operator 2D: primal";
    }
    else if (weak_form == "gauss points") {
      weak_sign = -1.0;
      pk_name = "PK operator 2D: gauss points";
    }
  } else {
    bool request_faces(true), request_edges(true);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx, request_faces, request_edges);
    pk_name = "PK operator 3D";
  }

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create global operator 
  // -- flux term
  ParameterList op_list = plist.sublist(pk_name).sublist("flux operator");
  op_list.sublist("schema").set<std::string>("dg basis", dg_basis);
  auto op_flux = Teuchos::rcp(new PDE_AdvectionRiemann(op_list, mesh));
  auto global_op = op_flux->global_operator();
  const WhetStone::DG_Modal& dg = op_flux->dg();

  // -- volumetric advection term
  op_list = plist.sublist(pk_name).sublist("advection operator");
  op_list.sublist("schema").set<std::string>("dg basis", dg_basis);
  auto op_adv = Teuchos::rcp(new PDE_Abstract(op_list, global_op));

  // -- reaction term
  op_list = plist.sublist(pk_name).sublist("reaction operator");
  op_list.sublist("schema").set<std::string>("dg basis", dg_basis);
  auto op_reac = Teuchos::rcp(new PDE_Reaction(op_list, global_op));
  double Kreac = op_list.get<double>("coef");

  // analytic solution
  int order = op_list.sublist("schema").get<int>("method order");
  int nk = (dim == 2) ? (order + 1) * (order + 2) / 2 :
                        (order + 1) * (order + 2) * (order + 3) / 6;

  AnalyticDG ana(mesh, order, true);

  // set up problem coefficients
  // -- reaction function
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);

  auto K = Teuchos::rcp(new CompositeVector(*cvs));
  auto Kc = K->ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    (*Kc)[0][c] = Kreac;
  }

  // -- velocity function
  auto velc = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost));
  auto velf = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));

  WhetStone::VectorPolynomial v;
  ana.VelocityTaylor(AmanziGeometry::Point(dim), 0.0, v); 
  
  for (int c = 0; c < ncells_wghost; ++c) {
    (*velc)[c] = v;
    (*velc)[c] *= -weak_sign;
  }

  for (int f = 0; f < nfaces_wghost; ++f) {
    (*velf)[f] = v * (mesh->face_normal(f) * weak_sign);
  }

  // -- divergence of velocity
  //    non-conservative formulation leads to Kn = Kreac - div(v)
  auto Kn = Teuchos::rcp(new std::vector<WhetStone::Polynomial>());
  WhetStone::Polynomial divv = Divergence(v);

  if (!conservative_form && weak_form == "dual") {
    auto tmp = divv;
    tmp *= -1.0;
    tmp(0, 0) += Kreac;

    Kn->resize(ncells_wghost, tmp);
  }

  // -- source term is calculated using method of manufactured solutions
  //    f = K p + div (v p) = K p + v . grad p + p div(v)
  WhetStone::Polynomial sol, src;
  WhetStone::Polynomial pc(dim, order);
  WhetStone::DenseVector data(pc.size());
  WhetStone::NumericalIntegration numi(mesh);

  Epetra_MultiVector& rhs_c = *global_op->rhs()->ViewComponent("cell");
  for (int c = 0; c < ncells; ++c) {
    const Point& xc = mesh->cell_centroid(c);

    v.ChangeOrigin(xc);
    divv.ChangeOrigin(xc);

    ana.SolutionTaylor(xc, 0.0, sol);

    src = Kreac * sol + v * Gradient(sol);
    if (conservative_form) src += divv * sol;

    for (auto it = pc.begin(); it < pc.end(); ++it) {
      int n = it.PolynomialPosition();

      WhetStone::Polynomial cmono(dim, it.multi_index(), 1.0);
      cmono.set_origin(xc);      

      WhetStone::Polynomial tmp = src * cmono;      

      data(n) = numi.IntegratePolynomialCell(c, tmp);
    }

    // -- convert moment to my basis
    dg.cell_basis(c).LinearFormNaturalToMy(data);
    for (int n = 0; n < pc.size(); ++n) {
      rhs_c[n][c] = data(n);
    }
  }

  // -- boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

  WhetStone::Polynomial coefs;

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    const Point& normal = mesh->face_normal(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[dim - 1]) < 1e-6 || fabs(xf[dim - 1] - 1.0) < 1e-6) {
      Point vp(v[0].Value(xf), v[1].Value(xf));

      if (vp * normal < -1e-12) {
        bc_model[f] = OPERATOR_BC_DIRICHLET;
        if (weak_sign < 0.0) bc_model[f] = OPERATOR_BC_DIRICHLET_TYPE2;

        ana.SolutionTaylor(xf, 0.0, coefs);
        data = coefs.coefs();

        for (int i = 0; i < nk; ++i) {
          bc_value[f][i] = data(i);
        }
      } else if (weak_sign < 0.0) {
        bc_model[f] = OPERATOR_BC_REMOVE;
      }
    }
  }

  // populate the global operator
  op_flux->SetBCs(bc, bc);
  op_flux->Setup(velc, velf);
  op_flux->UpdateMatrices(*velf);
  op_flux->ApplyBCs(true, true, true);

  op_adv->Setup(velc, false);
  op_adv->UpdateMatrices();

  if (conservative_form || weak_form == "primal" || weak_form == "gauss points")
    op_reac->Setup(Kc);
  else 
    op_reac->Setup(Kn, false);
  op_reac->UpdateMatrices(Teuchos::null, Teuchos::null);

  // create preconditoner
  global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"),
          "GMRES", plist.sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  CompositeVector& rhs = *global_op->rhs();
  CompositeVector solution(rhs);
  solution.PutScalar(0.0);

  global_op->ApplyInverse(rhs, solution);

  if (MyPID == 0) {
    std::cout << "dG solver (gmres): ||r||=" << global_op->residual() 
              << " itr=" << global_op->num_itrs()
              << " code=" << global_op->returned_code() 
              << " order=" << order << std::endl;

    // visualization
    const Epetra_MultiVector& p = *solution.ViewComponent("cell");
    GMV::open_data_file(*mesh, (std::string)"operators.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "solution");
    if (order > 0) {
      GMV::write_cell_data(p, 1, "gradx");
      GMV::write_cell_data(p, 2, "grady");
    }
    GMV::close_data_file();
  }

  CHECK(global_op->num_itrs() < 200);

  // compute solution error
  solution.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *solution.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int;
  ana.ComputeCellError(dg, p, 0.0, pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int);

  if (MyPID == 0) {
    sol.ChangeOrigin(AmanziGeometry::Point(dim));
    std::cout << "\nEXACT solution: " << sol << std::endl;
    printf("Mean:     L2(p)=%12.9f  Inf(p)=%12.9f  itr=%3d\n", pl2_mean, pinf_mean, global_op->num_itrs());
    printf("Total:    L2(p)=%12.9f  Inf(p)=%12.9f\n", pl2_err, pinf_err);
    printf("Integral: L2(p)=%12.9f\n", pl2_int);
    CHECK(pl2_err < 1e-10 && pinf_err < 1e-10);
  } 
}


TEST(OPERATOR_ADVECTION_STEADY_DG) {
  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "primal", false, "orthonormalized");
  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "primal", false, "normalized");
  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "primal", false, "regularized");
  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "dual", true, "normalized");
  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "dual", true, "regularized");
  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "dual", false);
  AdvectionSteady<AnalyticDG02>(3, "cubic", 3, "dual", true);

  AdvectionSteady<AnalyticDG03>(2, "test/median7x8.exo", 8, "gauss points", false, "orthonormalized");
}


