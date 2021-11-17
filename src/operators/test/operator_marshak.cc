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
#include "Mesh_MSTK.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionFactory.hh"
#include "UpwindFlux.hh"

#include "operator_marshak_testclass.hh"

void RunTestMarshak(std::string op_list_name, double TemperatureFloor) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->getRank();

  if (MyPID == 0) std::cout << "\nTest: Simulating nonlinear Marshak wave" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_marshak.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  // RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 3.0, 1.0, 200, 10);
  RCP<const Mesh> mesh = meshfactory.create("test/marshak.exo");

  // Create nonlinear coefficient.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh, TemperatureFloor));

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // modify diffusion coefficient
  CompositeVectorSpace K_map;
  K_map.SetMesh(mesh);
  K_map.AddComponent("cell", AmanziMesh::CELL, 1);
  auto K = Teuchos::rcp(new TensorVector(K_map));

  std::vector<WhetStone::Tensor<DefaultHostMemorySpace>> host_tensors(K->size());
  WhetStone::Tensor<DefaultHostMemorySpace> Kc;
  Kc.Init(2, 1);
  Kc(0,0) = 1.0;
  for (int c = 0; c < K->size(); c++) {
    host_tensors[c].assign(Kc);
  }
  K->Init(host_tensors); 
  
  // create boundary data (no mixed bc)
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  auto bc_model = bc->bc_model();
  auto bc_value = bc->bc_value();

  for (int f = 0; f < nfaces_wghost; f++) {
    const Point& xf = mesh->face_centroid(f);

    if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = 0.0;
    } else if (fabs(xf[0]) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = knc->TemperatureSource;
    } else if (fabs(xf[0] - 3.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = knc->TemperatureFloor;
    }
  }

  // Create and initialize solution (temperature) field.
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  if (op_list_name == "diffusion operator Sff") {
    cvs->SetMesh(mesh)->SetGhosted(true)
       ->AddComponent("cell", AmanziMesh::CELL, 1)
       ->AddComponent("face", AmanziMesh::FACE, 1);
  } else {
    cvs->SetMesh(mesh)->SetGhosted(true)
       ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(cvs->CreateSpace()));
  solution->putScalar(knc->TemperatureFloor);  // solution at time T=0

  // Create and initialize flux field.
 Teuchos::RCP<CompositeVectorSpace> cvs_flux = Teuchos::rcp(new CompositeVectorSpace());
  cvs_flux->SetMesh(mesh)
      ->SetComponent("face", AmanziMesh::FACE, 1)
      ->SetGhosted(true);
  auto flux = cvs_flux->Create();
  //op->UpdateFlux(solution.ptr(), flux.ptr());

  CompositeVector heat_capacity(cvs->CreateSpace());
  heat_capacity.putScalar(1.0);

  

  // Create upwind model
  ParameterList& ulist = plist.sublist("PK operator").sublist("upwind");
  UpwindFlux upwind(mesh);
  upwind.Init(ulist);

  // MAIN LOOP
  double tstop = plist.get<double>("simulation time", 0.5);
  int step(0);
  double snorm(0.0);
  
  double t(0.0), dt(1e-4);
  while (t < tstop) {
    solution->ScatterMasterToGhosted();


    // update bc
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->face_centroid(f);
      if (fabs(xf[0]) < 1e-6) bc_value[f] = knc->exact(t + dt, xf);
    }

    // upwind heat conduction coefficient
    knc->UpdateValues(*solution, *bc);
    upwind.Compute(*flux, *solution, bc_model, *knc->values());

    // add diffusion operator
    Teuchos::ParameterList olist = plist.sublist("PK operator").sublist(op_list_name);
    PDE_DiffusionFactory diff_factory;
    Teuchos::RCP<PDE_Diffusion> op = diff_factory.Create(olist, mesh, bc);
    op->Setup(K, knc->values(), knc->derivatives());
    op->UpdateMatrices(flux.ptr(), solution.ptr());
    // get the global operator
    Teuchos::RCP<Operator> global_op = op->global_operator();

    // add accumulation terms
    PDE_Accumulation op_acc(AmanziMesh::CELL, global_op);
    op_acc.AddAccumulationDelta(*solution, heat_capacity, heat_capacity, dt, "cell");

    // apply BCs and assemble
    op->ApplyBCs(true, true, true);
    global_op->set_inverse_parameters("Hypre AMG", plist.sublist("preconditioners"), "Amanzi GMRES", plist.sublist("solvers"));
    global_op->initializeInverse();
    global_op->computeInverse();

    auto sol_new = solution->ViewComponent("cell");
    auto sol_old(sol_new);


    CompositeVector rhs = *global_op->rhs();
    //    global_op->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);
    global_op->applyInverse(rhs, *solution);

    step++;
    t += dt;

    snorm = solution->norm2();

    if (MyPID == 0) {
      printf("%3d  ||r||=%11.6g  itr=%2d  ||sol||=%11.6g  t=%7.4f  dt=%7.4f\n",
          step, global_op->residual(), global_op->num_itrs(), snorm, t, dt);
    }

    // Change time step based on solution change.
    // We use empiric algorithm insired by Levenberg-Marquardt 
    auto sol_diff(sol_old);
    for(int i = 0 ; i < sol_diff.size(); ++i){
      sol_diff(i,0) = sol_diff(i,0)*-1.0+sol_new(i,0)*1.0; 
    }
    //sol_diff.update(1.0, sol_new, -1.0);

    double ds_rel(0.0);
    for (int c = 0; c < ncells_owned; c++) {
      ds_rel = std::max(ds_rel, sol_diff(c,0) / (1e-3 + sol_old(c,0) + sol_diff(c,0)));
    }
    double ds_rel_local = ds_rel;
    Teuchos::reduceAll(*(mesh->get_comm()),Teuchos::REDUCE_MAX,1,&ds_rel_local,&ds_rel); 

    if (ds_rel < 0.05) {
      dt *= 1.2;
    } else if (ds_rel > 0.10) {
      dt *= 0.8;
    }
    // dT = std::min(dT, 0.002);
  }

  // calculate errors
  const auto& p = solution->ViewComponent("cell");
  double pl2_err(0.0), pnorm(0.0);

  for (int c = 0; c < ncells_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double err = p(c,0) - knc->exact(t, xc);
    pl2_err += err * err;
    pnorm += p(c,0) * p(c,0);
  }
  double tmp = pl2_err;
  Teuchos::reduceAll(*(mesh->get_comm()),Teuchos::REDUCE_SUM,1,&tmp,
                     &pl2_err);

  //mesh->get_comm()->SumAll(&tmp, &pl2_err, 1);
  tmp = pnorm;
  Teuchos::reduceAll(*(mesh->get_comm()),Teuchos::REDUCE_SUM,1,&tmp,
                     &pnorm);
  //mesh->get_comm()->SumAll(&tmp, &pnorm, 1);

  pl2_err = std::pow(pl2_err / pnorm, 0.5);
  pnorm = std::pow(pnorm, 0.5);
  printf("||dp||=%10.6g  ||p||=%10.6g\n", pl2_err, pnorm);

  CHECK_CLOSE(0.0, pl2_err, 0.1);

  //if (MyPID == 0) {
  //  GMV::open_data_file(*mesh, (std::string)"operators.gmv");
  //  GMV::start_data();
  //  GMV::write_cell_data(p, 0, "solution");
  //  GMV::close_data_file();
  //}
}


/* *****************************************************************
* This test replaces tensor and boundary conditions by continuous
* functions. This is a prototype for heat conduction solvers.
* **************************************************************** */
// TEST(MARSHAK_NONLINEAR_WAVE_NLFV) {
//   RunTestMarshak("diffusion operator nlfv", 0.02);
// }

TEST(MARSHAK_NONLINEAR_WAVE_MFD) {
  RunTestMarshak("diffusion operator Sff", 0.0);
}

