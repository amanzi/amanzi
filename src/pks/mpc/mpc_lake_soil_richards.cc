#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"
#include "TreeOperator.hh"
#include "PDE_CouplingFlux.hh"
#include "InverseFactory.hh"
#include "PDE_Diffusion.hh"
#include "Operator_Schema.hh"

#include "mpc_lake_soil_richards.hh"

namespace Amanzi {

MPCLakeSoilRichards::MPCLakeSoilRichards(Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln) :
    PK(FElist, plist,  S, soln),
    StrongMPC<PK_BDF_Default>(FElist, plist,  S, soln) {}



void
MPCLakeSoilRichards::Setup(const Teuchos::Ptr<State>& S) {
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");

  std::cout << "In MPCLakeSoilRichards::Setup" << std::endl;

  // -- turn on coupling
  pks_list_->sublist(names[0]).set("coupled to soil via temp", true);
  pks_list_->sublist(names[1]).set("coupled to lake via temp", true);

  std::cout << "names[0] = " << names[0] << std::endl;
  std::cout << "names[1] = " << names[1] << std::endl;

//  // -- ensure local ops are suface ops
//  pks_list_->sublist(names[1]).sublist("diffusion preconditioner").set("surface operator", true);
//  pks_list_->sublist(names[1]).sublist("accumulation preconditioner").set("surface operator", true);

  domain_lake_ = pks_list_->sublist(names[0]).get<std::string>("domain name","lake");
  domain_soil_ = pks_list_->sublist(names[1]).get<std::string>("domain name","domain");

  std::cout << "domain_lake_ = " << domain_lake_ << std::endl;
  std::cout << "domain_soil_ = " << domain_soil_ << std::endl;

  // grab the meshes
  lake_mesh_ = S->GetMesh(domain_lake_);
  soil_mesh_ = S->GetMesh(domain_soil_);

  // cast the PKs
  lake_pk_ = sub_pks_[0];
  soil_pk_ = Teuchos::rcp_dynamic_cast<Amanzi::MPCCoupledSoil>(sub_pks_[1]);

  // require the coupling fields, claim ownership
  S->RequireField(Keys::getKey(domain_lake_,"lake_soil_temperature"), name_)
    ->SetMesh(lake_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // Create the preconditioner.
  // -- collect the preconditioners
//  precon_lake_ = lake_pk_->preconditioner();
//  precon_soil_ = soil_pk_->preconditioner();

  ti_list_ = Teuchos::sublist(plist_, "time integrator", true);

//  // -- set parameters for an inverse
//  Teuchos::ParameterList inv_list = plist_->sublist("inverse");
//  inv_list.setParameters(plist_->sublist("preconditioner"));
//  inv_list.setParameters(plist_->sublist("linear solver"));
//  precon_lake_->set_inverse_parameters(inv_list);

//  // -- push the surface local ops into the subsurface global operator
//  for (Operators::Operator::op_iterator op = precon_surf_->begin();
//       op != precon_surf_->end(); ++op) {
//    precon_lake_->OpPushBack(*op);
//  }

  // grab the debuggers
//  lake_db_ = lake_pk_->debugger();
//  soil_db_ = soil_pk_->debugger();
//  water_->set_db(domain_db_);

  // call the MPC's setup, which calls the sub-pk's setups
  StrongMPC<PK_BDF_Default>::Setup(S);

}

void MPCLakeSoilRichards::Initialize(const Teuchos::Ptr<State>& S) {

  auto comm = Amanzi::getDefaultComm();

   // initialize coupling terms
  S->GetFieldData(Keys::getKey(domain_lake_,"lake_soil_temperature"), name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_lake_,"lake_soil_temperature"), name_)->set_initialized();

  // Initialize all sub PKs.

  StrongMPC<PK_BDF_Default>::Initialize(S);

//  // // ensure continuity of ICs... surface takes precedence.
//  // CopySurfaceToSubsurface(*S->GetFieldData(Keys::getKey(domain_surf_,"pressure"), sub_pks_[1]->name()),
//  //       		  S->GetFieldData(Keys::getKey(domain_ss_,"pressure"), sub_pks_[0]->name()).ptr());
//  // ensure continuity of ICs... subsurface takes precedence.
//  CopySubsurfaceToSurface(*S->GetFieldData(Keys::getKey(domain_ss_,"pressure"), sub_pks_[0]->name()),
//			  S->GetFieldData(Keys::getKey(domain_surf_,"pressure"), sub_pks_[1]->name()).ptr());

//  // Initialize my timestepper.
//  PK_BDF_Default::Initialize(S);

  std::cout << "Solution_" << std::endl;

  solution_->Print(std::cout,false);

//  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
//  op_tree_lake_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

//  std::cout << "op_tree_lake_ = " << op_tree_lake_ << std::endl;
//
//  std::cout << "op_tree_lake_ Structure start:\n" << op_tree_lake_->PrintDiagnostics() << std::endl;

  auto op0 = lake_pk_->my_operator(Operators::OPERATOR_MATRIX)->Clone();
//  std::cout << "lake_pk_->my_operator(Operators::OPERATOR_MATRIX) = " << lake_pk_->my_operator(Operators::OPERATOR_MATRIX) << std::endl;
//  std::cout << "soil_pk_->my_operator(Operators::OPERATOR_MATRIX) = " << soil_pk_->my_operator(Operators::OPERATOR_MATRIX) << std::endl;
//  auto op1 = soil_pk_->my_operator(Operators::OPERATOR_MATRIX)->Clone();
//  std::cout << "soil_pk_->preconditioner() = " << soil_pk_->preconditioner() << std::endl;
//  auto op1 = soil_pk_->preconditioner()->Clone();
  auto op1 = soil_pk_->preconditioner()->Clone();

  std::cout << "beginning op0->A()" << std::endl;
  std::cout << "op0 = " << op0 << std::endl;
  op0->SymbolicAssembleMatrix();
  op0->AssembleMatrix();
  std::cout << *op0->A() << std::endl;

//  auto tvs = Teuchos::rcp(new TreeVectorSpace(op0->get_row_map()));
//
//  auto tvs_all = Teuchos::rcp(new TreeVectorSpace(comm));
//  tvs_all->PushBack(tvs);
//  tvs_all->PushBack(tvs);
//  tvs_all->PushBack(tvs);

//  op_tree_lake_ = Teuchos::rcp(new Operators::TreeOperator(tvs_all,tvs_all));

//  std::cout << "op_tree_lake_ = " << op_tree_lake_ << std::endl;

  std::cout << "Check 1 in mpc_lake_soil_richards.cc" << std::endl;

  auto op100 = soil_pk_->sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX);
  auto op111 = soil_pk_->sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX);

  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op0->DomainMap()))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(op1->RangeMap())));
//  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op100->DomainMap()))));
//  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(op111->DomainMap()))));
  op_tree_lake_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  std::cout << "op_tree_lake_ = " << op_tree_lake_ << std::endl;
  std::cout << "op_tree_lake_ Structure start:\n" << op_tree_lake_->PrintDiagnostics() << std::endl;

  op_tree_lake_->set_operator_block(0, 0, op0);
//  op_tree_lake_->set_operator_block(1, 1, op1);
  op_tree_lake_->set_block(1, 1, op1);

  std::cout << "op_tree_lake_ Structure with diagonal blocks:\n" << op_tree_lake_->PrintDiagnostics() << std::endl;

  std::cout << "Check 2 in mpc_lake_soil_richards.cc" << std::endl;

//  auto op100 = op1->get_operator_block(0,0);
//  auto op111 = op1->get_operator_block(1,1);

//  op1->set_operator_block(0,0,op100);
//  op1->set_operator_block(1,1,op111);

//  std::cout << "CHECK op_tree_lake_ :" << op_tree_lake_->get_operator() << std::endl;

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto mesh_lake = S_->GetMesh("lake");
  auto mesh_soil = S_->GetMesh("domain");

  std::cout << "Check 3 in mpc_lake_soil_richards.cc" << std::endl;

  std::cout << "solution_->SubVector(1)->SubVector(0)->Data() = " << solution_->SubVector(1)->SubVector(0)->Data() << std::endl;
  std::cout << "solution_->SubVector(1)->SubVector(1)->Data() = " << solution_->SubVector(1)->SubVector(1)->Data() << std::endl;

  auto& mmap10 = solution_->SubVector(1)->SubVector(0)->Data()->ViewComponent("cell", false)->Map();
  std::cout << "Check 40 in mpc_lake_soil_richards.cc" << std::endl;
  auto& gmap10 = solution_->SubVector(1)->SubVector(0)->Data()->ViewComponent("cell", true)->Map();
  auto& mmap11 = solution_->SubVector(1)->SubVector(1)->Data()->ViewComponent("cell", false)->Map();
  std::cout << "Check 41 in mpc_lake_soil_richards.cc" << std::endl;
  auto& gmap11 = solution_->SubVector(1)->SubVector(1)->Data()->ViewComponent("cell", true)->Map();
  std::cout << "Check 5 in mpc_lake_soil_richards.cc" << std::endl;
  int npoints_owned = mmap10.NumMyPoints()+mmap11.NumMyPoints();
  std::cout << "Check 6 in mpc_lake_soil_richards.cc" << std::endl;

  std::cout << "npoints_owned = " << npoints_owned << std::endl;

  auto cvs_lake = Teuchos::rcp(new CompositeVectorSpace());
  auto cvs_soil = Teuchos::rcp(new CompositeVectorSpace());

  cvs_lake->SetMesh(mesh_lake)->SetGhosted(true)
              ->AddComponent("cell", AmanziMesh::CELL, 1);

  cvs_soil->SetMesh(mesh_soil)->SetGhosted(true)
              ->AddComponent("cell", AmanziMesh::CELL, 1);


//  // -- indices transmissibimility coefficients for matrix-fracture flux
//  const auto& kn = *S_->GetFieldData("fracture-normal_permeability")->ViewComponent("cell");
//  double gravity;
//  S->GetConstantVectorData("gravity")->Norm2(&gravity);

  S_->GetFieldEvaluator("thermal_conductivity")->HasFieldChanged(S_.ptr(), "soil thermo");
  const auto& lambda_s = *S_->GetFieldData("thermal_conductivity")->ViewComponent("cell");

  S_->GetFieldEvaluator("lake-thermal_conductivity")->HasFieldChanged(S_.ptr(), "lake thermo");
  const auto& lambda_l = *S_->GetFieldData("lake-thermal_conductivity")->ViewComponent("cell");

//<---

//---> lake "0"

  auto& mmap1 = solution_->SubVector(0)->Data()->ViewComponent("cell", false)->Map();
  auto& gmap1 = solution_->SubVector(0)->Data()->ViewComponent("cell", true)->Map();
  int npoints_owned1 = mmap1.NumMyPoints();

  int ncells_owned_f1 = mesh_lake->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto inds_lake1 = std::make_shared<std::vector<std::vector<int> > >(npoints_owned1);
  auto inds_soil1 = std::make_shared<std::vector<std::vector<int> > >(npoints_owned1);
  auto values1 = std::make_shared<std::vector<double> >(npoints_owned1);

  auto values00 = std::make_shared<std::vector<double> >(npoints_owned1);

  auto inds_lake10 = std::make_shared<std::vector<std::vector<int> > >(1);
  auto inds_soil10 = std::make_shared<std::vector<std::vector<int> > >(1);
  auto values10 = std::make_shared<std::vector<double> >(1);

  std::cout << "npoints_owned1 = " << npoints_owned1 << std::endl;

  std::cout << "ncells_owned_f1 = " << ncells_owned_f1 << std::endl;

  int np1(0);
  for (int c = 0; c < ncells_owned_f1; ++c) {
    int f = mesh_lake->entity_get_parent(AmanziMesh::CELL, c);
    double area = mesh_lake->cell_volume(c);
    int first = mmap1.FirstPointInElement(f);
    int ndofs = mmap1.ElementSize(f);

    std::cout << "first = " << first << std::endl;
    std::cout << "ndofs = " << ndofs << std::endl;

    std::cout << "np1 = " << np1 << std::endl;

    for (int k = 0; k < ndofs; ++k) {
      (*inds_lake1)[np1].resize(1);
      (*inds_soil1)[np1].resize(1);

      (*values1)[np1] = 0.;
      (*values00)[np1] = 0.; //(lambda_l[0][c]+lambda_l[0][c]);
      if (c == 0) {
          (*values1)[np1] = -lambda_l[0][c]/area;
          (*values10)[0] = -lambda_l[0][c]/area;
          (*inds_lake10)[0].resize(1);
          (*inds_soil10)[0].resize(1);
          (*inds_soil10)[0][0] = ncells_owned_f1-1;
          (*inds_lake10)[0][0] = c;
      }   // normal??
      if (c == ncells_owned_f1-1) { (*values00)[np1] = 0.; //-lambda_l[0][c]/area; //-40.;
      (*inds_lake1)[np1][0] = c;
      (*inds_soil1)[np1][0] = 0;} //lambda_l[0][c]/area; }
      if (c == 0) { (*values00)[np1] = lambda_s[0][c]/area; //40.;
      (*inds_lake1)[np1][0] = c;
      (*inds_soil1)[np1][0] = 0;} //(lambda_l[0][c]+lambda_s[0][c])/area;}
//      if (c == 1) { (*values00)[np1] = -lambda_l[0][c]/area;}
//      std::cout << "lambda_l[0][c] = " << lambda_l[0][c] << std::endl;}
      np1++;
    }
  }

  inds_lake1->resize(np1);
  inds_soil1->resize(np1);
  values1->resize(np1);

  values00->resize(np1);

  inds_lake10->resize(1);
  inds_soil10->resize(1);
  values10->resize(1);

//<---

  //---> soil "1"

  int ncells_owned_f = mesh_soil->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto inds_lake = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  auto inds_soil = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  auto values = std::make_shared<std::vector<double> >(npoints_owned);

  auto values11 = std::make_shared<std::vector<double> >(npoints_owned);

  auto inds_lake01 = std::make_shared<std::vector<std::vector<int> > >(1);
  auto inds_soil01 = std::make_shared<std::vector<std::vector<int> > >(1);
  auto values01 = std::make_shared<std::vector<double> >(1);

  std::cout << "ncells_owned_f = " << ncells_owned_f << std::endl;

  int np(0);
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_soil->entity_get_parent(AmanziMesh::CELL, c);
    double area = mesh_soil->cell_volume(c);
    int first = mmap10.FirstPointInElement(f);
    int ndofs = mmap10.ElementSize(f);

    std::cout << "first = " << first << std::endl;
    std::cout << "ndofs = " << ndofs << std::endl;

    std::cout << "np = " << np << std::endl;

    for (int k = 0; k < ndofs; ++k) {
      (*inds_lake)[np].resize(1);
      (*inds_soil)[np].resize(1);

      (*values)[np] = 0.;
      (*values11)[np] = 0.;
      if (c == ncells_owned_f-1) {
          (*values)[np] = -lambda_s[0][c]/area;
          (*values01)[0] = -lambda_s[0][c]/area;
          (*inds_lake01)[0].resize(1);
          (*inds_soil01)[0].resize(1);
          (*inds_lake01)[0][0] = 0;
          (*inds_soil01)[0][0] = c;
      } // normal???
      if (c == 0) { (*values11)[np] = 0.; //-lambda_s[0][c]/area; //-40.;
            (*inds_lake)[np][0] = 0;
            (*inds_soil)[np][0] = c;}
      if (c == ncells_owned_f-1) { (*values11)[np] = lambda_l[0][c]/area; //40.;
            (*inds_lake)[np][0] = 0;
            (*inds_soil)[np][0] = c;} //lambda_s[0][c]/area;}//(lambda_l[0][c]+lambda_s[0][c])/area;}
//      if (c == ncells_owned_f-2) { (*values11)[np] = -lambda_s[0][c]/area;}
//      std::cout << "lambda_s[0][c] = " << lambda_s[0][c] << std::endl;}
      np++;
    }
  }

  inds_lake->resize(np);
  inds_soil->resize(np);
  values->resize(np);

  values11->resize(np);

  inds_lake01->resize(1);
  inds_soil01->resize(1);
  values01->resize(1);

  // -- operators
  Teuchos::ParameterList oplist;

  auto op_coupling00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_lake, cvs_lake, inds_lake1, inds_lake1, op0));
  op_coupling00->Setup(values00, 1.0);                            // values00 will be ADDED (with factor)
  op_coupling00->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto op_coupling01 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_lake, cvs_soil, inds_lake01, inds_soil01));
  op_coupling01->Setup(values01,1.0);
  op_coupling01->UpdateMatrices(Teuchos::null, Teuchos::null);

  std::cout << "op_coupling01 = " << op_coupling01 << std::endl;

  auto op_coupling10 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_soil, cvs_lake, inds_soil10, inds_lake10)); // cvs defines the block size
  op_coupling10->Setup(values10, 1.0);
  op_coupling10->UpdateMatrices(Teuchos::null, Teuchos::null);

  std::cout << "op_coupling10 = " << op_coupling10 << std::endl;


  auto op1_clone = op1->Clone();

  // DO COUPLING ON CLONE

  // ALTERNATIVELY:
  auto op11 = op1->get_operator_block(1,1)->Clone();

  auto op_coupling11 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_soil, cvs_soil, inds_soil, inds_soil, op11));
  op_coupling11->Setup(values11, 1.0);
  op_coupling11->UpdateMatrices(Teuchos::null, Teuchos::null);

  op1->set_operator_block(1, 1, op_coupling11->global_operator());

  op_tree_lake_->set_block(1, 1, op1);

  std::cout << "values00" << std::endl;
  for (int c = 0; c < (*values00).size(); c++) {
      std::cout << (*values00)[c] << std::endl;
  }

  std::cout << "values01" << std::endl;
  for (int c = 0; c < (*values01).size(); c++) {
      std::cout << (*values01)[c] << std::endl;
  }

  std::cout << "values10" << std::endl;
  for (int c = 0; c < (*values10).size(); c++) {
      std::cout << (*values10)[c] << std::endl;
  }

  std::cout << "values11" << std::endl;
  for (int c = 0; c < (*values11).size(); c++) {
      std::cout << (*values11)[c] << std::endl;
  }


  std::cout << "op_coupling01->global_operator() = " << op_coupling01->global_operator() << std::endl;
  std::cout << "op_coupling10->global_operator() = " << op_coupling10->global_operator() << std::endl;

//  auto tvs0 = Teuchos::rcp(new TreeVectorSpace(op0->get_row_map()));
//  auto tvs1 = Teuchos::rcp(new TreeVectorSpace(op0->get_row_map()));

//  auto op_soil = soil_pk_->sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX)->Clone();

  auto tvs0 = Teuchos::rcp(new TreeVectorSpace(op0->get_row_map()));
//  auto tvs1 = Teuchos::rcp(new TreeVectorSpace(op1->get_row_map()));
  auto tvs1 = Teuchos::rcp(new TreeVectorSpace(op1->RangeMap()));
 //  auto tvs1 = Teuchos::rcp(new TreeVectorSpace(solution_->SubVector(1)->Map()));

  auto tvs00 = Teuchos::rcp(new TreeVectorSpace(comm));
  tvs00->PushBack(tvs0);

  std::cout << "tvs00" << std::endl;
  tvs00->Print(std::cout);

  auto tvs01 = Teuchos::rcp(new TreeVectorSpace(comm));
  tvs01 = tvs1;
//  tvs01->PushBack(tvs1);
//  tvs01->PushBack(tvs1);

  std::cout << "tvs01" << std::endl;
  tvs01->Print(std::cout);

  auto tvs10 = Teuchos::rcp(new TreeVectorSpace(comm));
  tvs10 = tvs1;
//  tvs10->PushBack(tvs1);
//  tvs10->PushBack(tvs1);

  std::cout << "tvs10" << std::endl;
  tvs10->Print(std::cout);

  auto op01 = Teuchos::rcp(new Operators::TreeOperator(tvs00, tvs01));
  op01->set_operator_block(0, 1, op_coupling01->global_operator());
  op_tree_lake_->set_block(0, 1, op01);

  std::cout << "op01 Structure:\n" << op01->PrintDiagnostics() << std::endl;

  auto op10 = Teuchos::rcp(new Operators::TreeOperator(tvs01, tvs00)); // <--- should be tvs10, tvs00 but doesn't work
  op10->set_operator_block(1, 0, op_coupling10->global_operator());
  op_tree_lake_->set_block(1, 0, op10);

  std::cout << "op10 Structure:\n" << op10->PrintDiagnostics() << std::endl;

  // print capabilities
  std::cout << "Tree Operator Structure:\n" << op_tree_lake_->PrintDiagnostics() << std::endl;

  std::cout << "DONE test" << std::endl;
//  exit(0);

//  preconditioner_list_ = Teuchos::sublist(plist_, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(plist_, "solvers", true);
  ti_list_ = Teuchos::sublist(plist_, "time integrator", true);

//  std::string name = ti_list_->get<std::string>("preconditioner");
//  std::string ls_name = ti_list_->get<std::string>("preconditioner enhancement", "none");
//  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(name, *preconditioner_list_,
//                        ls_name, *linear_operator_list_,
//                        true);


  std::cout << "op0 Structure:\n" << op0->PrintDiagnostics() << std::endl;

  std::cout << "op1 Structure:\n" << op1->PrintDiagnostics() << std::endl;

  std::cout << "op0->A()" << std::endl;
  std::cout << "op0 = " << op0 << std::endl;
  op0->SymbolicAssembleMatrix();
  op0->AssembleMatrix();
  std::cout << *op0->A() << std::endl;

  std::cout << "op1->A()" << std::endl;
  std::cout << "op1 = " << op1 << std::endl;
  op1->SymbolicAssembleMatrix();
  op1->AssembleMatrix();
  std::cout << *op1->A() << std::endl;

  std::cout << "PRINT" << std::endl;
  std::cout << op_tree_lake_->PrintDiagnostics() << std::endl;
  std::cout << "PRINT DONE" << std::endl;

  op_tree_lake_->SymbolicAssembleMatrix();
  op_tree_lake_->AssembleMatrix();

  std::cout << "before A" << std::endl;
  std::cout << *op_tree_lake_->A() << std::endl;
  std::cout << "after A" << std::endl;
//  exit(0);


  std::cout << "plist_ = " << plist_ << std::endl;
  std::cout << "plist_->sublist('inverse') = " << plist_->sublist("inverse") << std::endl;
  op_tree_lake_->set_inverse_parameters(plist_->sublist("inverse"));
  op_tree_lake_->InitializeInverse();
//  std::cout << "inv_list = " << inv_list << std::endl;
//  op_tree_lake_->set_inverse_parameters(inv_list);
//  op_tree_lake_->InitializeInverse();

//  op_tree_lake_->ComputeInverse();

  // stationary solve is modelled with large dt. To pick the correct
  // boundary conditions, dt is negative. This assumes that we are at
  // the beginning of simulation.
//  if (ti_list_->isSublist("initialization")) {
//    // bool wells_on = ti_list_->sublist("initialization").get<bool>("active wells", false);
//    double dt(-1e+98), dt_solver;
//    bool fail = time_stepper_->TimeStep(dt, dt_solver, solution_);
//    if (fail) Exceptions::amanzi_throw("Solver for coupled Darcy flow did not converge.");
//  }

//  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
  if (true) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "MPC lake:" << std::endl
               << op_tree_lake_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl << std::endl;
  }


}

void
MPCLakeSoilRichards::set_states(const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<State>& S_inter,
                            const Teuchos::RCP<State>& S_next) {
  StrongMPC<PK_BDF_Default>::set_states(S,S_inter,S_next);
//  water_->set_states(S,S_inter,S_next);
}


// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk FunctionalResidual().
void
MPCLakeSoilRichards::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                            Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {

  std::cout << "FunctionalResidual START" << std::endl;

  // generate local matrices and apply sources and boundary conditions
  StrongMPC<PK_BDF_Default>::FunctionalResidual(t_old, t_new, u_old, u_new, f);



  // although, residual calculation can be completed using off-diagonal
  // blocks, we use global matrix-vector multiplication instead.
  op_tree_lake_->AssembleMatrix();
  int ierr = op_tree_lake_->ApplyAssembled(*u_new, *f);
  AMANZI_ASSERT(!ierr);

  // diagonal blocks in tree operator must be lake and soil
  AMANZI_ASSERT(sub_pks_[0]->name() == "lake thermo" &&
                  sub_pks_[1]->name() == "coupled soil");

  auto op0 = sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX);
//  auto op1 = sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX);
  auto op1 = soil_pk_->preconditioner();
//  auto op1 = soil_pk_->preconditioner()->Clone();

  std::cout << "FunctionalResidual op0->A()" << std::endl;
  std::cout << "op0 = " << op0 << std::endl;
  op0->SymbolicAssembleMatrix();
  op0->AssembleMatrix();
  std::cout << *op0->A() << std::endl;

  std::cout << "FunctionalResidual op1->A()" << std::endl;
  std::cout << "op1 = " << op1 << std::endl;
  op1->SymbolicAssembleMatrix();
  op1->AssembleMatrix();
  std::cout << *op1->A() << std::endl;
//
//  exit(0);

  soil_pk_ = Teuchos::rcp_dynamic_cast<Amanzi::MPCCoupledSoil>(sub_pks_[1]);
  auto op10 = soil_pk_->sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX);
  auto op11 = soil_pk_->sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX);

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pcA = soil_pk_->sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Operator> pcB = soil_pk_->sub_pks_[1]->preconditioner();

//  std::cout << "MPC LSR pcA" << std::endl;
//  pcA->SymbolicAssembleMatrix();
//  pcA->AssembleMatrix();
//  std::cout << *pcA->A() << std::endl;
//  std::cout << "MPC LSR pcB" << std::endl;
//  pcB->SymbolicAssembleMatrix();
//  pcB->AssembleMatrix();
//  std::cout << *pcB->A() << std::endl;

//  exit(0);


  f->SubVector(0)->Data()->Update(-1.0, *op0->rhs(), 1.0);
  f->SubVector(1)->SubVector(0)->Data()->Update(-1.0, *op10->rhs(), 1.0);
  f->SubVector(1)->SubVector(1)->Data()->Update(-1.0, *op11->rhs(), 1.0);
//  f->SubVector(1)->Data()->Update(-1.0, *op1->rhs(), 1.0);

//  sub_pks_[0]->FunctionalResidual(t_old, t_new, u_old->SubVector(0),
//            u_new->SubVector(0), f->SubVector(0));
//  soil_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
//            u_new->SubVector(1), f->SubVector(1));
//  soil_pk_->sub_pks_[0]->FunctionalResidual(t_old, t_new, u_old->SubVector(1)->SubVector(0),
//              u_new->SubVector(1)->SubVector(0), f->SubVector(1)->SubVector(0));
//  soil_pk_->sub_pks_[1]->FunctionalResidual(t_old, t_new, u_old->SubVector(1)->SubVector(1),
//                u_new->SubVector(1)->SubVector(1), f->SubVector(1)->SubVector(1));

  std::cout << *op_tree_lake_->A() << std::endl;

  u_old->Print(std::cout,false);

//  exit(0);

  // propagate updated info into state
//  Solution_to_State(*u_new, S_next_);
//
//  // Evaluate the surface flow residual
//  soil_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
//                            u_new->SubVector(1), g->SubVector(1));
//
////
////  // The residual of the surface flow equation provides the mass flux from
////  // subsurface to surface.
////
////  Epetra_MultiVector& source = *S_next_->GetFieldData(Keys::getKey(domain_surf_,"surface_subsurface_flux"),
////          name_)->ViewComponent("cell",false);
////  source = *g->SubVector(1)->Data()->ViewComponent("cell",false);
////
//  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
//  lake_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0),
//          u_new->SubVector(0), g->SubVector(0));
//
////  // All surface to subsurface fluxes have been taken by the subsurface.
////  g->SubVector(1)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  std::cout << "FunctionalResidual DONE" << std::endl;
}

// -- Apply preconditioner to u and returns the result in Pu.
int MPCLakeSoilRichards::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

  Pu->PutScalar(0.0);
  return op_tree_lake_->ApplyInverse(*u, *Pu);
//  int ierr = StrongMPC::ApplyPreconditioner(u,Pu);

  /*
  int ierr = 0;
  if (precon_type_ == PRECON_NONE) {
    *Pu = *u;
    ierr = 1;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    ierr = StrongMPC::ApplyPreconditioner(u,Pu);
  } else if (precon_type_ == PRECON_PICARD) {
    ierr = op_tree_lake_->ApplyInverse(*u, *Pu);
  } else if (precon_type_ == PRECON_EWC) {
    ierr = op_tree_lake_->ApplyInverse(*u, *Pu);
  }

  std::cout << "before ApplyInverse" << std::endl;
  ierr = op_tree_lake_->ApplyInverse(*u, *Pu);
  std::cout << "after ApplyInverse" << std::endl;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "PC * residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC*r_p"); vnames.push_back("  PC*r_T");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(Pu->SubVector(0)->Data().ptr());
    vecs.push_back(Pu->SubVector(1)->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }
  */


//  int ierr = StrongMPC::ApplyPreconditioner(u,Pu);

//  return (ierr > 0) ? 0 : 1;

}

// -- Update the preconditioner.
void
MPCLakeSoilRichards::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {

  std::cout << "UpdatePreconditioner START" << std::endl;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // order important -- subsurface's pk includes the surface's local ops, so
  // doing the subsurface 2nd re-inits the surface matrices (and doesn't
  // refill them).  This is why subsurface is first
  StrongMPC::UpdatePreconditioner(t, up, h);
//  op_tree_lake_->ComputeInverse();

  std::cout << "UpdatePreconditioner DONE" << std::endl;


}

// -- Modify the predictor.
bool
MPCLakeSoilRichards::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
//  bool modified = false;
//
//  // Calculate consistent faces
//  modified = soil_pk_->ModifyPredictor(h, u0->SubVector(0), u->SubVector(0));
//
////  // Merge surface cells with subsurface faces
////  if (modified) {
////
////    S_next_->GetFieldEvaluator(Keys::getKey(domain_surf_,"relative_permeability"))->HasFieldChanged(S_next_.ptr(),name_);
////    Teuchos::RCP<const CompositeVector> h_prev = S_inter_->GetFieldData(Keys::getKey(domain_surf_,"ponded_depth"));
////    MergeSubsurfaceAndSurfacePressure(*h_prev, u->SubVector(0)->Data().ptr(),
////            u->SubVector(1)->Data().ptr());
////  }
////
////
////  // Hack surface faces
////  bool newly_modified = false;
////  newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
////  newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
////  modified |= newly_modified;
////
////  // -- copy surf --> sub
////  if (newly_modified) {
////    CopySurfaceToSubsurface(*u->SubVector(1)->Data(), u->SubVector(0)->Data().ptr());
////  }
//
//  // Calculate consistent surface faces
//  modified |= lake_pk_->ModifyPredictor(h, u0->SubVector(1), u->SubVector(1));
//
//  return modified;

  bool modified(false);

  // potentially update faces
  modified = StrongMPC<PK_BDF_Default>::ModifyPredictor(h, u0, u);
  return modified;

}

// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCLakeSoilRichards::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();

//  // dump to screen
//  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
//    *vo_->os() << "NKA'd, PC'd correction." << std::endl;
//
//    std::vector<std::string> vnames;
//    vnames.push_back("p");
//    vnames.push_back("PC*p");
//
//    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
//    vecs.push_back(res->SubVector(0)->Data().ptr());
//    vecs.push_back(du->SubVector(0)->Data().ptr());
//
//    *vo_->os() << " Subsurface precon:" << std::endl;
//    domain_db_->WriteVectors(vnames, vecs, true);
//
//    vecs[0] = res->SubVector(1)->Data().ptr();
//    vecs[1] = du->SubVector(1)->Data().ptr();
//
//    *vo_->os() << " Surface precon:" << std::endl;
//    surf_db_->WriteVectors(vnames, vecs, true);
//  }

  // modify correction using sub-pk approaches
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified_res =
    StrongMPC<PK_BDF_Default>::ModifyCorrection(h, res, u, du);

//  // modify correction using water approaches
//  int n_modified = 0;
//  n_modified += water_->ModifyCorrection_WaterFaceLimiter(h, res, u, du);
//  double damping = water_->ModifyCorrection_WaterSpurtDamp(h, res, u, du);
//  n_modified += water_->ModifyCorrection_WaterSpurtCap(h, res, u, du, damping);

//  // -- accumulate globally
//  int n_modified_l = n_modified;
//  u->SubVector(0)->Data()->Comm()->SumAll(&n_modified_l, &n_modified, 1);
//  bool modified = (n_modified > 0) || (damping < 1.);

//  // -- calculate consistent subsurface cells
//  if (modified) {
//    // if (consistent_cells_) {
//    //   // Derive subsurface cell corrections.
//    //   precon_->UpdateConsistentCellCorrection(
//    //       *u->SubVector(0)->Data(),
//    //       du->SubVector(0)->Data().ptr());
//    // }
//
//    // Copy subsurface face corrections to surface cell corrections
//    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
//                            du->SubVector(1)->Data().ptr());
//  }

  // if (modified) {
  //   // Derive surface face corrections.
  //   UpdateConsistentFaceCorrectionWater_(res, du);
  // }

//  // dump to screen
//  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
//    *vo_->os() << "Modified correction." << std::endl;
//
//    std::vector<std::string> vnames;
//    vnames.push_back("p");
//    vnames.push_back("PC*p");
//
//    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
//    vecs.push_back(res->SubVector(0)->Data().ptr());
//    vecs.push_back(du->SubVector(0)->Data().ptr());
//
//    *vo_->os() << " Subsurface precon:" << std::endl;
//    domain_db_->WriteVectors(vnames, vecs, true);
//
//    vecs[0] = res->SubVector(1)->Data().ptr();
//    vecs[1] = du->SubVector(1)->Data().ptr();
//
//    *vo_->os() << " Surface precon:" << std::endl;
//    surf_db_->WriteVectors(vnames, vecs, true);
//  }

  return modified_res;

//  return (modified_res || modified) ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING :
//      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}


} // namespace
