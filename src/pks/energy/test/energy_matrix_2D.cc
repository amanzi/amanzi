/*
  Energy

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
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
#include "EnergyTwoPhase_PK.hh"
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "State.hh"
#include "WhetStoneDefs.hh"

/* **************************************************************** 
* Generates a preconditioner for the implicit discretization of
* the thermal operator.
* ************************************************************** */
TEST(ENERGY_2D_MATRIX) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::Energy;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: 2D homogeneous medium, preconditioner" << std::endl;

  // read parameter list 
  std::string xmlFileName = "test/energy_matrix_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  const Teuchos::RCP<Teuchos::ParameterList> plist = 
      Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  pref.push_back(Framework::STK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 17, 17);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  // initialize the Energy process kernel 
  Teuchos::ParameterList pk_tree = plist->sublist("PKs").sublist("energy");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  EnergyTwoPhase_PK* EPK = new EnergyTwoPhase_PK(pk_tree, plist, S, soln);
  EPK->Setup(S.ptr());
std::cout << "Passes EPK.Setup()" << std::endl;
  S->Setup();
std::cout << "Passed S.Setup()" << std::endl;
  S->InitializeFields();
std::cout << "Passed S.InitilizeFields()" << std::endl;
  S->InitializeEvaluators();
std::cout << "Passed S.InitilizeEvaluators()" << std::endl;
  EPK->Initialize(S.ptr());
std::cout << "Passed EPK.Initilize()" << std::endl;
  S->WriteDependencyGraph();
  S->CheckAllFieldsInitialized();

  // modify the default state for the problem at hand 
  // create the initial temperature function 
  std::string passwd("thermal");
  Epetra_MultiVector& temperature = *S->GetFieldData("temperature", passwd)->ViewComponent("cell");
  temperature.PutScalar(273.0);
  Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator("temperature"))->SetFieldAsChanged(S.ptr());

  // compute conductivity
  EPK->UpdateConductivityData(S.ptr());

  // create boundary data
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  Teuchos::RCP<BCs> bc = Teuchos::rcp(new BCs(mesh, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<double>& bc_value = bc->bc_value();
  
  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6) {
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = 200.0;
    }
  }
  
  // create diffusion operator 
  const Teuchos::ParameterList& elist = plist->sublist("PKs").sublist("energy");
  Teuchos::ParameterList oplist = elist.sublist("operators")
                                       .sublist("diffusion operator")
                                       .sublist("preconditioner");
  PDE_DiffusionFactory opfactory;
  Teuchos::RCP<PDE_Diffusion> op1 = opfactory.Create(oplist, mesh, bc);
  op1->SetBCs(bc, bc);

  // populate the diffusion operator
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(EPK->get_K());
  op1->Setup(Kptr, Teuchos::null, Teuchos::null);
  op1->UpdateMatrices(Teuchos::null, Teuchos::null);
  Teuchos::RCP<Operator> op = op1->global_operator();

  // add accumulation term
  Teuchos::RCP<PDE_Accumulation> op2 = Teuchos::rcp(new PDE_Accumulation(AmanziMesh::CELL, op));
  double dT = 1.0;
  CompositeVector solution(op->DomainMap());

  S->GetFieldEvaluator("energy")->HasFieldDerivativeChanged(S.ptr(), passwd, "temperature");
  const CompositeVector& dEdT = *S->GetFieldData("denergy_dtemperature");

  op2->AddAccumulationDelta(solution, dEdT, dEdT, dT, "cell");

  // add advection term: u = q_l n_l c_v
  // we do not upwind n_l c_v  in this test.
  S->GetFieldEvaluator("internal_energy_liquid")->HasFieldDerivativeChanged(S.ptr(), passwd, "temperature");
  const Epetra_MultiVector& c_v = *S->GetFieldData("dinternal_energy_liquid_dtemperature")
      ->ViewComponent("cell", true);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell", true);

  Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(op->DomainMap()));
  Epetra_MultiVector& q_l = *flux->ViewComponent("face");

  AmanziMesh::Entity_ID_List cells;

  AmanziGeometry::Point velocity(1e-4, 1e-4);
  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    q_l[0][f] = velocity * normal;
    
    mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    double tmp(0.0);
    for (int i = 0; i < ncells; i++) {
      int c = cells[i];
      tmp += n_l[0][c] * c_v[0][c];
    }
    q_l[0][f] *= tmp / ncells;
  }

  Teuchos::ParameterList alist;
  Teuchos::RCP<PDE_AdvectionUpwind> op3 = Teuchos::rcp(new PDE_AdvectionUpwind(alist, op));
  op3->SetBCs(bc, bc);
  op3->Setup(*flux);
  op3->UpdateMatrices(flux.ptr());

  // build the matrix
  op1->ApplyBCs(true, true, true);
  op3->ApplyBCs(true, true, true);
  op->SymbolicAssembleMatrix();
  op->AssembleMatrix();

  // make preconditioner
  // Teuchos::RCP<Operator> op3 = Teuchos::rcp(new Operator(*op2));

  Teuchos::ParameterList slist = plist->sublist("preconditioners").sublist("Hypre AMG");
  op->set_inverse_parameters(slist);
  op->InitializeInverse();
  op->ComputeInverse();

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(temperature, 0, "temperature");
    GMV::close_data_file();
  }

  delete EPK;
}
