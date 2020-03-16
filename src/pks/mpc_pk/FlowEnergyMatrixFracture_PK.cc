/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples flow and energy in matrix and fractures.
*/

#include "FlowEnergy_PK.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "primary_variable_field_evaluator.hh"
#include "TreeOperator.hh"

#include "FlowEnergyMatrixFracture_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
FlowEnergyMatrixFracture_PK::FlowEnergyMatrixFracture_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& glist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln)
  : glist_(glist), 
    Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ =  Teuchos::rcp(new VerboseObject("CoupledThermalFlow_PK", vlist));
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(glist, "PKs");
  if (pks_list->isSublist(name_)) {
    plist_ = Teuchos::sublist(pks_list, name_); 
  } else {
    std::stringstream messagestream;
    messagestream << "There is no sublist for PK " << name_ << "in PKs list\n";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(plist_, "time integrator", true);  
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void FlowEnergyMatrixFracture_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_domain_ = S->GetMesh();
  mesh_fracture_ = S->GetMesh("fracture");
  int dim = mesh_domain_->space_dimension();

  Teuchos::ParameterList& elist = S->FEList();

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs, so we need to define it here
  // -- pressure
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
  if (!S->HasField("pressure")) {
    *S->RequireField("pressure", "flow")->SetMesh(mesh_domain_)->SetGhosted(true) = *cvs;

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "pressure");
    auto eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator("pressure", eval);
  }

  if (!S->HasField("temperature")) {
    *S->RequireField("temperature", "thermal")->SetMesh(mesh_domain_)->SetGhosted(true) = *cvs;

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "temperature");
    auto eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator("temperature", eval);
  }

  // -- darcy flux
  if (!S->HasField("darcy_flux")) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S->RequireField("darcy_flux", "flow")->SetMesh(mesh_domain_)->SetGhosted(true) 
      ->SetComponent(name, AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux for fracture
  if (!S->HasField("fracture-darcy_flux")) {
    auto cvs2 = Operators::CreateNonManifoldCVS(mesh_fracture_);
    *S->RequireField("fracture-darcy_flux", "flow")->SetMesh(mesh_fracture_)->SetGhosted(true) = *cvs2;
  }

  // Require additional fields and evaluators
  Key normal_permeability_key_("fracture-normal_permeability");
  if (!S->HasField(normal_permeability_key_)) {
    S->RequireField(normal_permeability_key_, "state")->SetMesh(mesh_fracture_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // inform dependent PKs about coupling
  // -- flow 
  auto & mflow = glist_->sublist("PKs").sublist("flow matrix")
                        .sublist("physical models and assumptions");
  mflow.set<std::string>("coupled matrix fracture flow", "matrix");

  auto& fflow = glist_->sublist("PKs").sublist("flow fracture")
                       .sublist("physical models and assumptions");
  fflow.set<std::string>("coupled matrix fracture flow", "fracture");

  // -- energy
  auto& menergy = glist_->sublist("PKs").sublist("energy matrix")
                         .sublist("physical models and assumptions");
  menergy.set<std::string>("coupled matrix fracture energy", "matrix");

  auto& fenergy = glist_->sublist("PKs").sublist("energy fracture")
                         .sublist("physical models and assumptions");
  fenergy.set<std::string>("coupled matrix fracture energy", "fracture");

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup(S);
}


/* ******************************************************************* 
* Initialization create a tree operator to assemble global matrix
******************************************************************* */
void FlowEnergyMatrixFracture_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  PK_MPCStrong<PK_BDF>::Initialize(S);

  // diagonal blocks in tree operator are the FlowEnergy PKs
  auto pk_matrix = Teuchos::rcp_dynamic_cast<FlowEnergy_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<FlowEnergy_PK>(sub_pks_[1]);

  auto tvs = Teuchos::rcp(new TreeVectorSpace(solution_->Map()));
  op_tree_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  op_tree_->SetTreeOperatorBlock(0, 0, pk_matrix->op_tree());
  op_tree_->SetTreeOperatorBlock(2, 2, pk_fracture->op_tree());

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto mesh_matrix = S_->GetMesh("domain");
  auto mesh_fracture = S_->GetMesh("fracture");

  // auto& mmap = solution_->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  // auto& gmap = solution_->SubVector(0)->Data()->ViewComponent("face", true)->Map();
  // int npoints_owned = mmap.NumMyPoints();

  // auto cvs_matrix = Teuchos::rcp(new CompositeVectorSpace());
  // auto cvs_fracture = Teuchos::rcp(new CompositeVectorSpace());

  // cvs_matrix->SetMesh(mesh_matrix)->SetGhosted(true)
  //           ->AddComponent("face", AmanziMesh::FACE, Teuchos::rcpFromRef(mmap), Teuchos::rcpFromRef(gmap), 1);

  // cvs_fracture->SetMesh(mesh_matrix)->SetGhosted(true)
  //             ->AddComponent("cell", AmanziMesh::CELL, 1);

  // // -- indices transmissibimility coefficients for matrix-fracture flux
  // const auto& kn = *S_->GetFieldData("fracture-normal_permeability")->ViewComponent("cell");
  // double rho = *S->GetScalarData("fluid_density");
  // double gravity;
  // S->GetConstantVectorData("gravity")->Norm2(&gravity);

  // int ncells_owned_f = mesh_fracture->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  // auto inds_matrix = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  // auto inds_fracture = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  // auto values = std::make_shared<std::vector<double> >(npoints_owned);

  // int np(0);
  // for (int c = 0; c < ncells_owned_f; ++c) {
  //   int f = mesh_fracture->entity_get_parent(AmanziMesh::CELL, c);
  //   double area = mesh_fracture->cell_volume(c);
  //   int first = mmap.FirstPointInElement(f);
  //   int ndofs = mmap.ElementSize(f);

  //   for (int k = 0; k < ndofs; ++k) {
  //     (*inds_matrix)[np].resize(1);
  //     (*inds_fracture)[np].resize(1);
  //     (*inds_matrix)[np][0] = first + k;
  //     (*inds_fracture)[np][0] = c;

  //     (*values)[np] = kn[0][c] * area / gravity;
  //     np++;
  //   }
  // }

  // inds_matrix->resize(np);
  // inds_fracture->resize(np);
  // values->resize(np);

  // // -- operators
  // Teuchos::ParameterList oplist;

  // auto op_coupling00 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
  //     oplist, cvs_matrix, cvs_matrix, inds_matrix, inds_matrix, pk_matrix->op()));
  // op_coupling00->Setup(values, 1.0);
  // op_coupling00->UpdateMatrices(Teuchos::null, Teuchos::null);

  // auto op_coupling01 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
  //     oplist, cvs_matrix, cvs_fracture, inds_matrix, inds_fracture));
  // op_coupling01->Setup(values, -1.0);
  // op_coupling01->UpdateMatrices(Teuchos::null, Teuchos::null);

  // auto op_coupling10 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
  //     oplist, cvs_fracture, cvs_matrix, inds_fracture, inds_matrix));
  // op_coupling10->Setup(values, -1.0);
  // op_coupling10->UpdateMatrices(Teuchos::null, Teuchos::null);

  // auto op_coupling11 = Teuchos::rcp(new Operators::PDE_CouplingFlux(
  //     oplist, cvs_fracture, cvs_fracture, inds_fracture, inds_fracture, pk_fracture->op()));
  // op_coupling11->Setup(values, 1.0);
  // op_coupling11->UpdateMatrices(Teuchos::null, Teuchos::null);  

  // op_tree_->SetOperatorBlock(0, 1, op_coupling01->global_operator());
  // op_tree_->SetOperatorBlock(1, 0, op_coupling10->global_operator());

  // // create a global problem
  // pk_matrix->op_diff()->ApplyBCs(true, true, true);

  op_tree_->SymbolicAssembleMatrix();
  op_tree_->AssembleMatrix();

  // Test SPD properties of the matrix.
  // VerificationTV ver(op_tree_);
  // ver.CheckMatrixSPD();

  // stationary solve is modelled with large dt. To pick the correct
  // boundary conditions, dt is negative. This assumes that we are at
  // the beginning of simulation.
  bool fail;
  double dt(-1e+98), dt_solver;
  fail = time_stepper_->TimeStep(dt, dt_solver, solution_);

  if (fail) Exceptions::amanzi_throw("Solver for coupled Darcy flow did not converge.");
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowEnergyMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  return fail;
}


/* ******************************************************************* 
* Residual evaluation
******************************************************************* */
void FlowEnergyMatrixFracture_PK::FunctionalResidual(
    double t_old, double t_new,
    Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
    Teuchos::RCP<TreeVector> f)
{
  // generate local matrices and apply sources and boundary conditions
  // although, residual calculation can be completed using off-diagonal
  // blocks, we use global matrix-vector multiplication instead.
  PK_MPCStrong<PK_BDF>::FunctionalResidual(t_old, t_new, u_old, u_new, f);
  op_tree_->AssembleMatrix();

  int ierr = op_tree_->ApplyAssembled(*u_new, *f);
  AMANZI_ASSERT(!ierr);
  
  // // diagonal blocks in tree operator and the Darcy PKs
  // auto pk_matrix = Teuchos::rcp_dynamic_cast<Flow::Darcy_PK>(sub_pks_[0]);
  // auto pk_fracture = Teuchos::rcp_dynamic_cast<Flow::Darcy_PK>(sub_pks_[1]);

  // f->SubVector(0)->Data()->Update(-1, *pk_matrix->op()->rhs(), 1);
  // f->SubVector(1)->Data()->Update(-1, *pk_fracture->op()->rhs(), 1);
}


/* ******************************************************************* 
* Preconditioner update
******************************************************************* */
void FlowEnergyMatrixFracture_PK::UpdatePreconditioner(
    double t, Teuchos::RCP<const TreeVector> up, double h)
{
  std::string name = ti_list_->get<std::string>("preconditioner");
  Teuchos::ParameterList pc_list = preconditioner_list_->sublist(name);

  op_tree_->InitPreconditioner(pc_list);
}


/* ******************************************************************* 
* Application of preconditioner
******************************************************************* */
int FlowEnergyMatrixFracture_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                                     Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_tree_->ApplyInverse(*X, *Y);
}

}  // namespace Amanzi

