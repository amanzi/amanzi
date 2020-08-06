/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples Transport in matrix and fracture
  using implicit scheme.
*/

#include "LinearOperatorGMRES.hh"
#include "Op_Diagonal.hh"
#include "PK_BDF.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "TransportImplicit_PK.hh"
#include "TreeOperator.hh"
#include "UniqueLocalIndex.hh"

#include "TransportMatrixFractureImplicit_PK.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
TransportMatrixFractureImplicit_PK::TransportMatrixFractureImplicit_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& glist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln)
  : glist_(glist),
    Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  auto pk_list = Teuchos::sublist(glist, "PKs", true);
  tp_list_ = Teuchos::sublist(pk_list, pk_name, true);

  Teuchos::ParameterList vlist;
  vo_ = Teuchos::rcp(new VerboseObject("TranCoupledImplicit_PK", vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void TransportMatrixFractureImplicit_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_domain_ = S->GetMesh();
  mesh_fracture_ = S->GetMesh("fracture");

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs

  // -- darcy flux in matrix
  if (!S->HasField("darcy_flux")) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S->RequireField("darcy_flux", "state")->SetMesh(mesh_domain_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux in fracture
  if (!S->HasField("fracture-darcy_flux")) {
    auto cvs = Operators::CreateNonManifoldCVS(mesh_fracture_);
    *S->RequireField("fracture-darcy_flux", "state")->SetMesh(mesh_fracture_)->SetGhosted(true) = *cvs;
  }

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup(S);
}


/* ******************************************************************* 
* Initialization creates a tree operator to assemble global matrix
******************************************************************* */
void TransportMatrixFractureImplicit_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  PK_MPCStrong<PK_BDF>::Initialize(S);

  // set a huge time step that will be limited by advance step
  set_dt(1e+98);

  TimestepControllerFactory<TreeVector> factory;
  auto ts_list = tp_list_->sublist("time integrator").sublist("BDF1");
  ts_control_ = factory.Create(ts_list, Teuchos::null, Teuchos::null);

  // diagonal blocks in tree operator are the Transport Implicit PKs
  auto pk_matrix = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[1]);

  auto tvs = Teuchos::rcp(new TreeVectorSpace(my_solution_->Map()));
  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  op_tree_matrix_->SetOperatorBlock(0, 0, pk_matrix->op());
  op_tree_matrix_->SetOperatorBlock(1, 1, pk_fracture->op());

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto cvs_matrix = Teuchos::rcp(new CompositeVectorSpace());
  auto cvs_fracture = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix->SetMesh(mesh_domain_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);
  cvs_fracture->SetMesh(mesh_fracture_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- indices are fluxes on matrix-fracture interface
  int ncells_owned_f = mesh_fracture_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto inds_matrix = std::make_shared<std::vector<std::vector<int> > >(2 * ncells_owned_f);
  auto inds_fracture = std::make_shared<std::vector<std::vector<int> > >(2 * ncells_owned_f);
  auto values = std::make_shared<std::vector<double> >(2 * ncells_owned_f, 0.0);

  int np(0);
  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->entity_get_parent(AmanziMesh::CELL, c);
    mesh_domain_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    AMANZI_ASSERT(ncells == 2);

    for (int k = 0; k < ncells; ++k) {
      (*inds_matrix)[np].resize(1);
      (*inds_fracture)[np].resize(1);
      (*inds_matrix)[np][0] = cells[k];
      (*inds_fracture)[np][0] = c;
      np++;
    }
  }

  // -- operators
  Teuchos::ParameterList oplist;

  op_coupling00_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_matrix, cvs_matrix, inds_matrix, inds_matrix, pk_matrix->op()));
  op_coupling00_->Setup(values, 1.0);
  op_coupling00_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling01_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_matrix, cvs_fracture, inds_matrix, inds_fracture));
  op_coupling01_->Setup(values, -1.0);
  op_coupling01_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling10_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_fracture, cvs_matrix, inds_fracture, inds_matrix));
  op_coupling10_->Setup(values, -1.0);
  op_coupling10_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling11_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, cvs_fracture, cvs_fracture, inds_fracture, inds_fracture, pk_fracture->op()));
  op_coupling11_->Setup(values, 1.0);
  op_coupling11_->UpdateMatrices(Teuchos::null, Teuchos::null);  

  op_tree_matrix_->SetOperatorBlock(0, 1, op_coupling01_->global_operator());
  op_tree_matrix_->SetOperatorBlock(1, 0, op_coupling10_->global_operator());

  // create a global problem
  pk_matrix->op_adv()->ApplyBCs(true, true, true);
  pk_fracture->op_adv()->ApplyBCs(true, true, true);

  op_tree_matrix_->SymbolicAssembleMatrix();

  // Test SPD properties of the matrix.
  // VerificationTV ver(op_tree_);
  // ver.CheckMatrixSPD();
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK is complete." 
               << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool TransportMatrixFractureImplicit_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;

  auto pk_matrix = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[0]);
  auto pk_fracture = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[1]);

  // update coupling terms
  int ncells_owned_f = mesh_fracture_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto values1 = std::make_shared<std::vector<double> >(2 * ncells_owned_f, 0.0);
  auto values2 = std::make_shared<std::vector<double> >(2 * ncells_owned_f, 0.0);

  int np(0), dir, shift;
  AmanziMesh::Entity_ID_List cells;
  const auto& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face");
  const auto& mmap = flux.Map();

  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->entity_get_parent(AmanziMesh::CELL, c);
    int first = mmap.FirstPointInElement(f);

    mesh_domain_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    mesh_domain_->face_normal(f, false, cells[0], &dir);
    shift = Operators::UniqueIndexFaceToCells(*mesh_domain_, f, cells[0]);

    for (int k = 0; k < ncells; ++k) {
      // since cells are ordered differenty than points, we need a map
      double tmp = flux[0][first + shift] * dir;

      if (tmp > 0) 
        (*values1)[np] = tmp;
      else 
        (*values2)[np] = -tmp;

      dir = -dir;
      shift = 1 - shift;
      np++;
    }
  }

  // update accumulation term and the right-hand side
  const auto& phi_m = *S_->GetFieldData("porosity");
  const auto& phi_f = *S_->GetFieldData("fracture-porosity");
  const auto& tcc_m = *S_->GetFieldData("total_component_concentration");
  const auto& tcc_f = *S_->GetFieldData("fracture-total_component_concentration");

  // we assume that all components are aquesous
  int num_aqueous = tcc_m.ViewComponent("cell")->NumVectors();

  for (int i = 0; i < num_aqueous; i++) {
    auto tv_one = ExtractComponent_(tcc_m, tcc_f, i);

    pk_matrix->UpdateBoundaryData(t_old, t_new, i);
    pk_fracture->UpdateBoundaryData(t_old, t_new, i);

    pk_matrix->op()->rhs()->PutScalar(0.0);
    pk_matrix->op_acc()->local_op(0)->Rescale(0.0);
    pk_matrix->op_acc()->AddAccumulationDelta(*tv_one->SubVector(0)->Data(), phi_m, phi_m, dt, "cell");

    pk_fracture->op()->rhs()->PutScalar(0.0);
    pk_fracture->op_acc()->local_op(0)->Rescale(0.0);
    pk_fracture->op_acc()->AddAccumulationDelta(*tv_one->SubVector(1)->Data(), phi_f, phi_f, dt, "cell");

    // assemble the operators
    pk_matrix->op_adv()->UpdateMatrices(S_->GetFieldData("darcy_flux").ptr());
    pk_matrix->op_adv()->ApplyBCs(true, true, true);

    pk_fracture->op_adv()->UpdateMatrices(S_->GetFieldData("fracture-darcy_flux").ptr());
    pk_fracture->op_adv()->ApplyBCs(true, true, true);

    op_coupling00_->Setup(values1, 1.0);
    op_coupling00_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling01_->Setup(values2, -1.0);
    op_coupling01_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling10_->Setup(values1, -1.0);
    op_coupling10_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling11_->Setup(values2, 1.0);
    op_coupling11_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_tree_matrix_->AssembleMatrix();

    // create preconditioner
    auto pc_list = glist_->sublist("preconditioners").sublist("Hypre AMG");
    op_tree_matrix_->InitPreconditioner(pc_list);

    // create solver
    Teuchos::ParameterList slist = glist_->sublist("solvers")
                                          .sublist("GMRES with Hypre AMG").sublist("gmres parameters");
    AmanziSolvers::LinearOperatorGMRES<Operators::TreeOperator, TreeVector, TreeVectorSpace>
        solver(op_tree_matrix_, op_tree_matrix_);
    solver.Init(slist);
    solver.add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);

    TreeVector rhs(*tv_one), tv_aux(*tv_one);
    *rhs.SubVector(0)->Data() = *pk_matrix->op()->rhs();
    *rhs.SubVector(1)->Data() = *pk_fracture->op()->rhs();

    int ierr = solver.ApplyInverse(rhs, tv_aux);

    SaveComponent_(tv_aux, my_solution_, i);

    // process error code
    bool fail = (ierr == 0);
    if (fail) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Step failed." << std::endl;
      return fail;
    }
  }

  // output 
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    pk_matrix->VV_PrintSoluteExtrema(*my_solution_->SubVector(0)->Data()->ViewComponent("cell"), dt, " (m)");
    pk_fracture->VV_PrintSoluteExtrema(*my_solution_->SubVector(1)->Data()->ViewComponent("cell"), dt, " (f)");
  }

  dt_ = ts_control_->get_timestep(dt_, 1);

  return false;
}


/* ******************************************************************
* Spezialized implementation of implicit transport
****************************************************************** */
void TransportMatrixFractureImplicit_PK::CommitStep(
    double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  *S->GetFieldData("total_component_concentration", "state") = *my_solution_->SubVector(0)->Data();
  *S->GetFieldData("fracture-total_component_concentration", "state") = *my_solution_->SubVector(1)->Data();
}


/* ******************************************************************
* Create a copy of the i-th component
****************************************************************** */
Teuchos::RCP<TreeVector> TransportMatrixFractureImplicit_PK::ExtractComponent_(
    const CompositeVector& tcc_m, const CompositeVector& tcc_f, int component)
{
  auto tv = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> tv1 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> tv2 = Teuchos::rcp(new TreeVector());
  tv->PushBack(tv1);
  tv->PushBack(tv2);

  CompositeVectorSpace cvs_m, cvs_f;
  cvs_m.SetMesh(mesh_domain_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);
  cvs_f.SetMesh(mesh_fracture_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  auto cv0 = Teuchos::rcp(new CompositeVector(cvs_m));
  tv->SubVector(0)->SetData(cv0);
  *(*cv0->ViewComponent("cell"))(0) = *(*tcc_m.ViewComponent("cell"))(component);

  auto cv1 = Teuchos::rcp(new CompositeVector(cvs_f));
  tv->SubVector(1)->SetData(cv1);
  *(*cv1->ViewComponent("cell"))(0) = *(*tcc_f.ViewComponent("cell"))(component);

  return tv;
}


/* ******************************************************************
* Save the i-th component to a solution vector
****************************************************************** */
void TransportMatrixFractureImplicit_PK::SaveComponent_(
    const TreeVector& tv_one, const Teuchos::RCP<TreeVector>& tv_all, int component)
{
  for (int i = 0; i < 2; ++i) {
    *(*tv_all->SubVector(i)->Data()->ViewComponent("cell"))(component) = 
        *(*tv_one.SubVector(i)->Data()->ViewComponent("cell"))(0);
  }
}

}  // namespace Amanzi

