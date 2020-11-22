/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Advection on non-manifold fracture network.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IterativeMethodGMRES.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "OutputXDMF.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwindDFN.hh"
#include "UniqueLocalIndex.hh"


/* *****************************************************************
* TBW.
* **************************************************************** */
void RunTest(double gravity) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: Transport in fracture network, gravity=" << gravity << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_diffusion_dfn.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an SIMPLE mesh framework
  ParameterList region_list = plist->sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  auto mlist = Teuchos::sublist(plist, "mesh", true);
  MeshFactory meshfactory(comm, gm, mlist);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create("test/fractures.exo");

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // create Darcy flux
  auto cvsf = Operators::CreateNonManifoldCVS(mesh);
  auto flux = Teuchos::rcp(new CompositeVector(*cvsf));
  Epetra_MultiVector& flux_f = *flux->ViewComponent("face", true);
  const auto& map = flux->Map().Map("face", true);

  int dir;
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point v(1.0, 0.0, 1.0);
  for (int f = 0; f < nfaces_owned; ++f) {
    int g = map->FirstPointInElement(f);
    int ndofs = map->ElementSize(f);

    mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (ndofs > 1) CHECK(ndofs == cells.size());

    for (int i = 0; i < ndofs; ++i) {
      int c = cells[i];
      auto normal = mesh->face_normal(f, false, c, &dir); 
      normal *= dir;  // natural normal
      
      int g2 = g + Operators::UniqueIndexFaceToCells(*mesh, f, c);
      flux_f[0][g2] = v * normal;
    }
  }

  // create boundary data
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = 1.0;
    }
  }

  // create solution 
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

  CompositeVector solution(*cvs), solution_new(*cvs);
  solution.PutScalar(0.0);

  // create advection operator
  Teuchos::ParameterList olist = plist->sublist("PK operator").sublist("advection operator");
  Teuchos::RCP<Operators::PDE_AdvectionUpwindDFN> op_adv = Teuchos::rcp(new PDE_AdvectionUpwindDFN(olist, mesh));
  Teuchos::RCP<Operator> global_op = op_adv->global_operator();

  // add accumulation operator 
  double dt(0.1);
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, global_op));
  op_acc->AddAccumulationDelta(solution, dt, "cell");

  // populate advection operator
  op_adv->Setup(*flux);
  op_adv->SetBCs(bc, bc);
  op_adv->UpdateMatrices(flux.ptr());

  // apply BCs and assemble
  op_adv->ApplyBCs(true, true, true);
    
  // create inverse
  global_op->set_inverse_parameters("Hypre AMG", plist->sublist("preconditioners"),
                               "GMRES", plist->sublist("solvers"));
  global_op->InitializeInverse();
  global_op->ComputeInverse();

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  // time stepping
  double t(0.0);
  for (int nstep = 0; nstep < 5; ++nstep) {

    CompositeVector& rhs = *global_op->rhs();
    global_op->ApplyInverse(rhs, solution_new);

    // -- modify right-hand side and solution
    const auto& old_c = *solution.ViewComponent("cell");
    const auto& new_c = *solution_new.ViewComponent("cell");
    const auto& rhs_c = *rhs.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; ++c) 
      rhs_c[0][c] += (new_c[0][c] - old_c[0][c]) * mesh->cell_volume(c) / dt;

    solution = solution_new;

    double a;
    rhs.Norm2(&a);

    io.InitializeCycle(t, nstep, "");
    io.WriteVector(*new_c(0), "solution", AmanziMesh::CELL);
    io.FinalizeCycle();

    // verify solution bounds and monotone decrease away from sources at x=0
    for (int c = 0; c < ncells_owned; ++c) {
      const auto& xc = mesh->cell_centroid(c);
      CHECK(new_c[0][c] <= 1.0);
      for (int c2 = 0; c2 < ncells_owned; ++c2) {
        const auto& xc2 = mesh->cell_centroid(c2);
        if (xc2[0] - xc[0] > 1e-6 && fabs(xc2[1] - xc[1]) < 1e-6
                                  && fabs(xc2[2] - xc[2]) < 1e-6) CHECK(new_c[0][c2] <= new_c[0][c]); 
      }
    }
  }
}


TEST(TRANSPORT_IN_FRACTURES) {
  RunTest(0.0);
}


