/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Tests for FCT.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "LeastSquare.hh"
#include "MeshFactory.hh"

// Amanzi::Operators
#include "FCT.hh"
#include "LimiterCell.hh"
#include "OperatorDefs.hh"
#include "ReconstructionCellGrad.hh"

double fun_field(const Amanzi::AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 return x*x*y + 2*x*y*y*y;
}

Amanzi::AmanziGeometry::Point fun_velocity(const Amanzi::AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 return Amanzi::AmanziGeometry::Point(2*x, -2*y);  // divergence-free
}


/* *****************************************************************
* Flux corrected transport in 2D
***************************************************************** */
std::pair<double, double> RunTest(int n) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  // int MyPID = comm->MyPID();
  // if (MyPID == 0) std::cout << "\nTest: Accuracy of FCT in 2D." << std::endl;

  // create rectangular mesh
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, n, n);

  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // initialize data
  // -- field
  auto field = Teuchos::rcp(new Epetra_MultiVector(mesh->cell_map(true), 1));

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    (*field)[0][c] = fun_field(xc);
  }

  // -- flux
  int dir;
  AmanziMesh::Entity_ID_List cells;

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)->SetGhosted(false)->AddComponent("face", AmanziMesh::FACE, 1);

  auto velocity = Teuchos::rcp(new CompositeVector(cvs));
  auto& velocity_f = *velocity->ViewComponent("face");

  auto flux_lo = Teuchos::rcp(new CompositeVector(cvs));
  auto flux_ho = Teuchos::rcp(new CompositeVector(cvs));

  auto flux_exact = Teuchos::rcp(new CompositeVector(cvs));
  auto flux_numer = Teuchos::rcp(new CompositeVector(cvs));

  auto& flux_lo_f = *flux_lo->ViewComponent("face");
  auto& flux_ho_f = *flux_ho->ViewComponent("face");
  auto& flux_exact_f = *flux_exact->ViewComponent("face");

  double dt(0.03 / n);
  for (int f = 0; f < nfaces_owned; f++) {
    mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    if (cells.size() == 2) {
      const auto& xf = mesh->face_centroid(f);
      const auto& normal = mesh->face_normal(f, false, cells[0], &dir);
      const auto& xc = (dir == 1) ? mesh->cell_centroid(cells[0])
                                  : mesh->cell_centroid(cells[1]);

      // low-order and high-order fluxes are from 1st to 2nd cell
      double tmp = (fun_velocity(xf) * normal) * dt;
      velocity_f[0][f] = tmp * dir;
 
      flux_lo_f[0][f] = tmp * fun_field(xc);  // default upwind
      flux_ho_f[0][f] = tmp * fun_field(xf);

      flux_exact_f[0][f] = flux_ho_f[0][f];
    }
  }

  // -- boundary conditions
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  auto& bc_model = bc->bc_model();
  auto& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(1.0 - xf[0]) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(1.0 - xf[1]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = fun_field(xf);
    }
  }

  // compute reconstruction
  Teuchos::ParameterList plist;
  plist.set<int>("polynomial_order", 1)
       .set<bool>("limiter extension for transport", false)
       .set<std::string>("limiter", "Barth-Jespersen")
       .set<std::string>("limiter stencil", "cell to all cells");

  auto lifting = Teuchos::rcp(new ReconstructionCellGrad(mesh));
  lifting->Init(plist);
  lifting->Compute(field);

  // create limiter (Do we need it?)
  auto limiter = Teuchos::rcp(new LimiterCell(mesh));
  limiter->Init(plist, velocity->ViewComponent("face"));

  // FCT
  FCT fct(mesh, mesh, limiter, field);
  fct.Compute(*flux_lo, *flux_ho, *bc, *flux_numer);

  // Error calculations
  double err[1];
  std::pair<double, double> out;

  flux_lo->Update(1.0, *flux_exact, -1.0);
  flux_lo->Norm2(&err[0]);
  out.first = err[0] * n;

  flux_numer->Update(1.0, *flux_exact, -1.0);
  flux_numer->Norm2(&err[0]);
  out.second = err[0] * n;

  std::cout << "error lower order: " << out.first << std::endl;
  std::cout << "error high order: " << out.second << std::endl;
  std::cout << "mean limiter: " << fct.get_alpha_mean() << std::endl;

  return out;
}


TEST(FCT_2D) {
  auto a1 = RunTest(10);
  auto a2 = RunTest(20);
  auto a3 = RunTest(40);

  std::vector<double> h({ 1.0/10, 1.0/30, 1.0/40 });
  std::vector<double> err1({ a1.first, a2.first, a3.first });
  double rate = Amanzi::Utils::bestLSfit(h, err1);
  CHECK(rate < 1.0); 
  std::cout << "\nError convergence rate: " << rate << std::endl;

  double err2 = std::max({ a1.second, a2.second, a3.second });
  CHECK(err2 < 2.0e-4); 
  std::cout << "Deviation from high-order flux: " << err2 << std::endl;
}

