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
#include "MatrixMarket_Tpetra.hpp"

// Amanzi
//#include "LinearOperatorGMRES.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "Tensor.hh"

// Amanzi::Operators
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionFactory.hh"

#include "operator_marshak_testclass.hh"

#define UPDATE_TEST 0

void writeMarshakMatrix(std::string op_list_name, double floor, bool jac) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int rank = comm->getRank();
  int size = comm->getSize();

  if (rank == 0) std::cout << "\nTest: Simulating nonlinear Marshak wave" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_marshak.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  
  // create an MSTK mesh framework
  ParameterList region_list = plist.sublist("regions");

  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, region_list, *comm));
  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 1);

  // Create nonlinear problem.
  Teuchos::RCP<HeatConduction> knc = Teuchos::rcp(new HeatConduction(mesh, floor));
  double t = 2.0;
  
  // set the diffusion coefficient
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nbfaces = mesh->num_entities(AmanziMesh::BOUNDARY_FACE, AmanziMesh::Parallel_type::OWNED);

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
  auto bc = Teuchos::rcp(new Operators::BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  {
    auto bc_value = bc->bc_value();
    auto bc_model = bc->bc_model();
    
    const auto& bf_map = *mesh->map(AmanziMesh::BOUNDARY_FACE, false);
    const auto& f_map = *mesh->map(AmanziMesh::FACE, false);
      
    for (int bf=0; bf!=bf_map.getLocalNumElements(); ++bf) {
      auto f = f_map.getLocalElement(bf_map.getGlobalElement(bf));
      auto xf = mesh->face_centroid(f);

      if (fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = 0.0;
      } else if (fabs(xf[0]) < 1e-6) {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = knc->exact(t, xf);
      } else if (fabs(xf[0] - 1.0) < 1e-6) {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value[f] = knc->exact(t, xf);
      } else {
        bc_model[f] = Operators::OPERATOR_BC_NONE;
      }
    }
  }

  // Create and initialize solution (temperature) field.
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  if (op_list_name != "fv: default") {
    cvs->AddComponent("face", AmanziMesh::FACE, 1);
  }

  auto solution = cvs->Create();
  {
    auto soln_c = solution->ViewComponent<Amanzi::MirrorHost>("cell", false);
    for (int c = 0; c < ncells_owned; ++c) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid_host(c);
      soln_c(c,0) = knc->exact(t, xc);
    }
  }
  if (solution->HasComponent("face")) {
    auto soln_f = solution->ViewComponent<Amanzi::MirrorHost>("face", false);
    for (int f = 0; f < nfaces_owned; ++f) {
      const AmanziGeometry::Point& xf = mesh->face_centroid_host(f);
      soln_f(f,0) = knc->exact(t, xf);
    }
  }
  solution->ScatterMasterToGhosted();
  solution->print(std::cout);

  // scalar coefficient
  knc->UpdateValues(*solution, *bc);

  // upwind
  knc->values()->ScatterMasterToGhosted("cell");
  knc->derivatives()->ScatterMasterToGhosted("cell");
  solution->ScatterMasterToGhosted("cell");
  
  {
    auto u_c = solution->ViewComponent<Amanzi::MirrorHost>("cell", true);
    auto kc = knc->values()->ViewComponent<Amanzi::MirrorHost>("cell", true);
    auto dkc = knc->derivatives()->ViewComponent<Amanzi::MirrorHost>("cell", true);

    auto kf = knc->values()->ViewComponent<Amanzi::MirrorHost>("face", false);
    auto kbf = knc->values()->ViewComponent<Amanzi::MirrorHost>("dirichlet_faces", false);

    auto dkf = knc->derivatives()->ViewComponent<Amanzi::MirrorHost>("face", false);
    auto dkbf = knc->derivatives()->ViewComponent<Amanzi::MirrorHost>("dirichlet_faces", false);

    auto bc_value = bc->bc_value();
    auto bc_model = bc->bc_model();
    
    for (int f=0; f!=nfaces_owned; ++f) {
      AmanziMesh::Entity_ID_View cells;
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
      if (cells.size() == 2) {
        kf(f,0) = u_c(cells(0),0) > u_c(cells(1),0) ? kc(cells(0),0) : kc(cells(1),0);
        dkf(f,0) = u_c(cells(0),0) > u_c(cells(1),0) ? dkc(cells(0),0) : dkc(cells(1),0);
      } else {
        if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
          int bf = mesh->exterior_face_map(false)->getLocalElement(mesh->face_map(false)->getGlobalElement(f));
          kf(f,0) = bc_value[f] > u_c(cells(0),0) ? kbf(bf,0) : kc(cells(0),0);
          dkf(f,0) = bc_value[f] > u_c(cells(0),0) ? dkbf(bf,0) : dkc(cells(0),0);
        } else {
          kf(f,0) = kc(cells(0),0);
          dkf(f,0) = dkc(cells(0),0);
        }
      }
    }
  }    
  // knc->values()->Print(std::cout);
  // knc->derivatives()->Print(std::cout);

  
  // create diffusion operator
  Teuchos::ParameterList olist;
  olist.set("discretization primary", op_list_name);
  olist.set("gravity", false);
  olist.set("nonlinear coefficient", "upwind: face");
  if (jac) {
    if (op_list_name == "fv: default") {
      olist.set("Newton correction", "true Jacobian");
    } else {
      olist.set("Newton correction", "approximate Jacobian");
    }
  }
  PDE_DiffusionFactory diff_factory;
  auto op = diff_factory.Create(olist, mesh, bc);
  op->Setup(K, knc->values(), knc->derivatives());
  op->UpdateMatrices(Teuchos::null, solution.ptr());

  // flux vector
  auto cvs_flux = Teuchos::rcp(new CompositeVectorSpace());
  cvs_flux->SetMesh(mesh)
      ->SetComponent("face", AmanziMesh::FACE, 1)
      ->SetGhosted(true);
  auto flux = cvs_flux->Create();
  op->UpdateFlux(solution.ptr(), flux.ptr());
  if (jac) op->UpdateMatricesNewtonCorrection(flux.ptr(), solution.ptr());

  // get the global operator
  Teuchos::RCP<Operator> global_op = op->global_operator();

  // apply BCs and assemble
  op->ApplyBCs(true, true, true);
  global_op->SymbolicAssembleMatrix();
  global_op->AssembleMatrix();

  std::replace(op_list_name.begin(), op_list_name.end(), ' ', '_');
  std::replace(op_list_name.begin(), op_list_name.end(), ':', '_');

  
  std::stringstream fname_stream;
  fname_stream << "test/operator_marshak_assembly_";
  if (jac) fname_stream << "jac_";
  fname_stream << op_list_name << "_np" << size;

  std::string filename_test = fname_stream.str() + ".test";
  std::string filename_gold = fname_stream.str() + ".gold";

#if UPDATE_TEST
  filename_test = filename_gold;
#endif


  // write to disk
  global_op->WriteMatrix(filename_test);

  if (rank == 0) {
    // only test on rank 0, this gets ugly...
    Tpetra::MatrixMarket::Reader<Matrix_type> reader;
    auto comm_self = Amanzi::getCommSelf();

    auto test_mat = reader.readSparseFile((filename_test+".dat").c_str(), comm_self);
    auto gold_mat = reader.readSparseFile(filename_gold.c_str(), comm_self);

    auto test_map = test_mat->getRowMap();
    auto gold_map = gold_mat->getRowMap();

    for (int i=0; i!=test_map->getLocalNumElements(); ++i) {
      int test_n_cols = -1;
      Kokkos::View<const double*, Layout, Amanzi::HostSpace, Kokkos::MemoryManaged> test_vals; 
      Kokkos::View<const int*, Layout, Amanzi::HostSpace, Kokkos::MemoryManaged> test_cols;
      test_mat->getLocalRowView(i,test_cols,test_vals);
      test_n_cols = test_cols.size(); 
      
      int gold_n_cols = -1;
      Kokkos::View<const double*, Layout, Amanzi::HostSpace, Kokkos::MemoryManaged> gold_vals; 
      Kokkos::View<const int*, Layout, Amanzi::HostSpace, Kokkos::MemoryManaged> gold_cols;
      gold_mat->getLocalRowView(i, gold_cols, gold_vals);
      gold_n_cols = gold_cols.size(); 
      
      std::cout << i << ": ";
      CHECK_EQUAL(gold_n_cols, test_n_cols);
      if (test_n_cols == gold_n_cols) {
        for (int j=0; j!=gold_n_cols; ++j) {
          int k=0;
          for (; k!=test_n_cols; ++k) {
            if (gold_cols[j] == test_cols[k]) break;
          }
          CHECK(k < test_n_cols);
          if (k == test_n_cols) continue;
          
          std::cout << "(" << gold_cols[j] << ": " << gold_vals[j] << " == " << test_vals[k] << ")";
          CHECK_CLOSE(gold_vals[j], test_vals[j], 1.e-8);
        }
      }
      std::cout << std::endl;
    }
  }
}


/* *****************************************************************
* This test replaces tensor and boundary conditions by continuous
* functions. This is a prototype for heat conduction solvers.
* **************************************************************** */
TEST(MARSHAK_NONLINEAR_WAVE_FV) {
  writeMarshakMatrix("fv: default", 0.0, false);
}
TEST(MARSHAK_NONLINEAR_WAVE_FV_JAC) {
  writeMarshakMatrix("fv: default", 0.0, true);
}


TEST(MARSHAK_NONLINEAR_WAVE_MFD_TPFA) {
  writeMarshakMatrix("mfd: two-point flux approximation", 0.0, false);
}
TEST(MARSHAK_NONLINEAR_WAVE_MFD_TPFA_JAC) {
  writeMarshakMatrix("mfd: two-point flux approximation", 0.0, true);
}

// TEST(MARSHAK_NONLINEAR_WAVE_MFD) {
//   writeMarshakMatrix("mfd: default", 0.0, false);
// }
// TEST(MARSHAK_NONLINEAR_WAVE_MFD_JAC) {
//   writeMarshakMatrix("mfd: default", 0.0, true);
// }


