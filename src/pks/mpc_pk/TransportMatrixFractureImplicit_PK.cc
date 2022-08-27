/*
  This is the mpc_pk component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  Process kernel that couples Transport in matrix and fracture
  using implicit scheme.
*/

#include "InverseFactory.hh"
#include "Op_Diagonal.hh"
#include "PDE_CouplingFlux.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "TransportImplicit_PK.hh"
#include "TreeOperator.hh"
#include "UniqueLocalIndex.hh"

#include "TransportMatrixFractureImplicit_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
TransportMatrixFractureImplicit_PK::TransportMatrixFractureImplicit_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& glist,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln)
   : Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
     Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln),
     glist_(glist),
     soln_(soln)
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
void TransportMatrixFractureImplicit_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  matrix_vol_flowrate_key_ = "volumetric_flow_rate";
  fracture_vol_flowrate_key_ = "fracture-volumetric_flow_rate";

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs

  // -- darcy flux in matrix
  if (!S_->HasRecord(matrix_vol_flowrate_key_)) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>(matrix_vol_flowrate_key_, Tags::DEFAULT, "state")
      .SetMesh(mesh_domain_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux in fracture
  if (!S_->HasRecord(fracture_vol_flowrate_key_)) {
    auto cvs = Operators::CreateManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>(fracture_vol_flowrate_key_, Tags::DEFAULT, "state")
      .SetMesh(mesh_fracture_)->SetGhosted(true) = *cvs;
  }

  S_->Require<CV_t, CVS_t>("fracture-normal_diffusion", Tags::DEFAULT, "state")
    .SetMesh(mesh_fracture_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup();
}


/* *******************************************************************
* Initialization creates a tree operator to assemble global matrix
******************************************************************* */
void TransportMatrixFractureImplicit_PK::Initialize()
{
  PK_MPCStrong<PK_BDF>::Initialize();

  // set a huge time step that will be limited by advance step
  set_dt(1e+98);

  TimestepControllerFactory<TreeVector> factory;
  auto ts_list = tp_list_->sublist("time integrator").sublist("BDF1");
  ts_control_ = factory.Create(ts_list, Teuchos::null, Teuchos::null);

  // diagonal blocks in tree operator are the Transport Implicit PKs
  pk_matrix_ = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[0]);
  pk_fracture_ = Teuchos::rcp_dynamic_cast<Transport::TransportImplicit_PK>(sub_pks_[1]);
  AMANZI_ASSERT(pk_matrix_->domain() == "domain");

  auto tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(pk_matrix_->op()->get_row_map())));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(pk_fracture_->op()->get_row_map())));
  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  auto op0 = pk_matrix_->op()->Clone();
  auto op1 = pk_fracture_->op()->Clone();
  op_tree_matrix_->set_operator_block(0, 0, op0);
  op_tree_matrix_->set_operator_block(1, 1, op1);

  // off-diagonal blocks represent coupling terms
  fia_ = Teuchos::rcp(new FractureInsertion(mesh_domain_, mesh_fracture_));
  fia_->InitMatrixCellToFractureCell();

  // -- advection
  Teuchos::ParameterList oplist;

  op_coupling00_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, fia_->get_cvs_matrix(), fia_->get_cvs_matrix(), 
              fia_->get_inds_matrix(), fia_->get_inds_matrix(), op0));
  op_coupling00_->Setup(fia_->get_values(), 1.0);
  op_coupling00_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling01_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, fia_->get_cvs_matrix(), fia_->get_cvs_fracture(),
              fia_->get_inds_matrix(), fia_->get_inds_fracture()));
  op_coupling01_->Setup(fia_->get_values(), -1.0);
  op_coupling01_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling10_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, fia_->get_cvs_fracture(), fia_->get_cvs_matrix(),
              fia_->get_inds_fracture(), fia_->get_inds_matrix()));
  op_coupling10_->Setup(fia_->get_values(), -1.0);
  op_coupling10_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_coupling11_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
      oplist, fia_->get_cvs_fracture(), fia_->get_cvs_fracture(),
              fia_->get_inds_fracture(), fia_->get_inds_fracture(), op1));
  op_coupling11_->Setup(fia_->get_values(), 1.0);
  op_coupling11_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_tree_matrix_->set_operator_block(0, 1, op_coupling01_->global_operator());
  op_tree_matrix_->set_operator_block(1, 0, op_coupling10_->global_operator());

  // -- dispersion/diffusion
  flag_dispersion_ = false;
  for (const auto& pk : sub_pks_) {
    flag_dispersion_ |= Teuchos::rcp_dynamic_cast<Transport::Transport_PK>(pk)->get_flag_dispersion();
  }

  if (flag_dispersion_) {
    auto& cvs0 = pk_matrix_->op()->get_domain_map();
    auto mmap = cvs0->Map("face", false);
    auto gmap = cvs0->Map("face", true);

    fid_ = Teuchos::rcp(new FractureInsertion(mesh_domain_, mesh_fracture_)); 
    fid_->InitMatrixFaceToFractureCell(mmap, gmap);

    op_coupling00d_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
        oplist, fid_->get_cvs_matrix(), fid_->get_cvs_matrix(), 
                fid_->get_inds_matrix(), fid_->get_inds_matrix(), op0));
    op_coupling00d_->Setup(fid_->get_values(), 1.0);
    op_coupling00d_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling01d_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
        oplist, fid_->get_cvs_matrix(), fid_->get_cvs_fracture(),
                fid_->get_inds_matrix(), fid_->get_inds_fracture(), op_coupling01_->global_operator()));
    op_coupling01d_->Setup(fid_->get_values(), -1.0);
    op_coupling01d_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling10d_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
        oplist, fid_->get_cvs_fracture(), fid_->get_cvs_matrix(),
                fid_->get_inds_fracture(), fid_->get_inds_matrix(), op_coupling10_->global_operator()));
    op_coupling10d_->Setup(fid_->get_values(), -1.0);
    op_coupling10d_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling11d_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(
        oplist, fid_->get_cvs_fracture(), fid_->get_cvs_fracture(),
                fid_->get_inds_fracture(), fid_->get_inds_fracture(), op1));
    op_coupling11d_->Setup(fid_->get_values(), 1.0);
    op_coupling11d_->UpdateMatrices(Teuchos::null, Teuchos::null);
  }

  // create a global problem
  pk_matrix_->op_adv()->ApplyBCs(true, true, true);
  pk_fracture_->op_adv()->ApplyBCs(true, true, true);

  Teuchos::ParameterList& ti_list = tp_list_->sublist("time integrator");
  std::string name = ti_list.get<std::string>("preconditioner", "Hypre AMG");
  std::string ls_name = ti_list.get<std::string>("linear solver", "GMRES with Hypre AMG");
  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
      name, glist_->sublist("preconditioners"), ls_name, glist_->sublist("solvers"), true);
  inv_list.setName(name);
  op_tree_matrix_->set_inverse_parameters(inv_list);
  op_tree_matrix_->InitializeInverse();

  if (!flag_dispersion_)
    InitializeCVField(S_, *vo_, "fracture-normal_diffusion", Tags::DEFAULT, "state", 0.0);

  // time integrators
  int nspace_m, nspace_f, ntime;
  pk_matrix_->get_discretization_order(&nspace_m, &ntime);
  pk_fracture_->get_discretization_order(&nspace_f, &ntime);

  if (nspace_m == 2 && nspace_f == 2) {
    Teuchos::ParameterList& bdf1_list = ti_list.sublist("BDF1");
    bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  }

  // Test SPD properties of the matrix.
  // VerificationTV ver(op_tree_);
  // ver.CheckMatrixSPD();
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "matrix:" << std::endl
               << op_tree_matrix_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete."
               << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************* 
* Maybe we need a separate PK instead of the if-clause
******************************************************************* */
bool TransportMatrixFractureImplicit_PK::AdvanceStep(
    double t_old, double t_new, bool reinit) 
{
  bool fail;

  // make copy of primary unknowns, i.e. solution 
  auto& tcc_m = S_->GetW<CV_t>("total_component_concentration", Tags::DEFAULT, "state");
  auto& tcc_f = S_->GetW<CV_t>("fracture-total_component_concentration", Tags::DEFAULT, "state");

  num_aqueous_ = tcc_m.ViewComponent("cell")->NumVectors();
  CompositeVector tcc_m_copy(tcc_m), tcc_f_copy(tcc_f);

  fia_->SetValues(S_->Get<CV_t>(matrix_vol_flowrate_key_));

  // fork between low-order and high-order
  int nspace_m, nspace_f, ntime, tot_itrs;
  pk_matrix_->get_discretization_order(&nspace_m, &ntime);
  pk_fracture_->get_discretization_order(&nspace_f, &ntime);

  if (nspace_m == 1 && nspace_f == 1) {
    fail = AdvanceStepLO_(t_old, t_new, &tot_itrs);
  } else if (nspace_m == 2 && nspace_f == 2) {
    fail = AdvanceStepHO_(t_old, t_new, &tot_itrs);
  } else {
    AMANZI_ASSERT(false);
  }

  if (fail) {
    tcc_m = tcc_m_copy;
    tcc_f = tcc_f_copy;
    // *my_solution_ = my_solution_copy;
    ChangedSolution();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed: recover (" << my_solution_->size() << ") primary field" << std::endl;
  }

  // output
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    double dt = t_new - t_old;
    Teuchos::OSTab tab = vo_->getOSTab();
    pk_matrix_->VV_PrintSoluteExtrema(*tcc_m.ViewComponent("cell"), dt, " (m)");
    pk_fracture_->VV_PrintSoluteExtrema(*tcc_f.ViewComponent("cell"), dt, " (f)");
  }

  return fail;
}


/* *******************************************************************
* One time step for aqueous components only
******************************************************************* */
bool TransportMatrixFractureImplicit_PK::AdvanceStepLO_(
    double t_old, double t_new, int* tot_itrs)
{
  auto& tcc_m = S_->GetW<CV_t>("total_component_concentration", Tags::DEFAULT, "state");
  auto& tcc_f = S_->GetW<CV_t>("fracture-total_component_concentration", Tags::DEFAULT, "state");

  for (int i = 0; i < num_aqueous_; i++) {
    pk_matrix_->UpdateLinearSystem(t_old, t_new, i);
    pk_fracture_->UpdateLinearSystem(t_old, t_new, i);

    op_coupling00_->Setup(fia_->get_values(), 1.0);
    op_coupling00_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling01_->Setup(fia_->get_values2(), -1.0);
    op_coupling01_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling10_->Setup(fia_->get_values(), -1.0);
    op_coupling10_->UpdateMatrices(Teuchos::null, Teuchos::null);

    op_coupling11_->Setup(fia_->get_values2(), 1.0);
    op_coupling11_->UpdateMatrices(Teuchos::null, Teuchos::null);

    // assemble dispersion/diffusion operators
    if (flag_dispersion_) {
      const auto& kn = *S_->Get<CV_t>("fracture-normal_diffusion").ViewComponent("cell");
      fid_->SetValues(kn, 1.0);

      op_coupling00d_->Setup(fid_->get_values(), 1.0);
      op_coupling00d_->UpdateMatrices(Teuchos::null, Teuchos::null);

      op_coupling01d_->Setup(fid_->get_values(), -1.0);
      op_coupling01d_->UpdateMatrices(Teuchos::null, Teuchos::null);

      op_coupling10d_->Setup(fid_->get_values(), -1.0);
      op_coupling10d_->UpdateMatrices(Teuchos::null, Teuchos::null);

      op_coupling11d_->Setup(fid_->get_values(), 1.0);
      op_coupling11d_->UpdateMatrices(Teuchos::null, Teuchos::null);
    }

    // create solver
    op_tree_matrix_->AssembleMatrix();
    op_tree_matrix_->ComputeInverse();
 
    auto& tvs = op_tree_matrix_->DomainMap();
    TreeVector rhs_one(tvs), sol_one(tvs);
    *rhs_one.SubVector(0)->Data() = *pk_matrix_->op()->rhs();
    *rhs_one.SubVector(1)->Data() = *pk_fracture_->op()->rhs();

    int ierr = op_tree_matrix_->ApplyInverse(rhs_one, sol_one);

    *(*tcc_m.ViewComponent("cell"))(i) = *(*sol_one.SubVector(0)->Data()->ViewComponent("cell"))(0);
    *(*tcc_f.ViewComponent("cell"))(i) = *(*sol_one.SubVector(1)->Data()->ViewComponent("cell"))(0);

    bool fail = (ierr != 0);
    if (fail) {
      dt_ = ts_control_->get_timestep(dt_, -1);
      return fail;
    }
  }

  dt_ = ts_control_->get_timestep(dt_, 1);
  return false;
}


/* ******************************************************************* 
*
******************************************************************* */
bool TransportMatrixFractureImplicit_PK::AdvanceStepHO_(
    double t_old, double t_new, int* tot_itrs) 
{
  double dt_next;

  dt_ = t_new - t_old;
  *tot_itrs = bdf1_dae_->number_nonlinear_steps();

  auto tcc_m = S_->GetPtrW<CV_t>("total_component_concentration", Tags::DEFAULT, "state");
  auto tcc_f = S_->GetPtrW<CV_t>("fracture-total_component_concentration", Tags::DEFAULT, "state");

  num_aqueous_ = tcc_m->ViewComponent("cell")->NumVectors();
  AMANZI_ASSERT(num_aqueous_ == 1);

  for (int i = 0; i < num_aqueous_; i++) {
    pk_matrix_->set_current_component(i);
    pk_fracture_->set_current_component(i);

    *soln_->SubVector(0)->Data() = *tcc_m;
    *soln_->SubVector(1)->Data() = *tcc_f;

    bool fail = bdf1_dae_->TimeStep(dt_, dt_next, soln_);
    dt_ = dt_next;
    if (fail) return fail;
  }

  *tcc_m = *soln_->SubVector(0)->Data();
  *tcc_f = *soln_->SubVector(1)->Data();

  bdf1_dae_->CommitSolution(dt_, soln_);
  *tot_itrs = bdf1_dae_->number_nonlinear_steps() - *tot_itrs;

  return false;
}


/* ******************************************************************
* Data were copied in the advance step
****************************************************************** */
void TransportMatrixFractureImplicit_PK::CommitStep(
    double t_old, double t_new, const Tag& tag)
{
}


/* ******************************************************************* 
* Modify defualt behavior:
******************************************************************* */
void TransportMatrixFractureImplicit_PK::CalculateDiagnostics(const Tag& tag)
{  
  int nspace_m, nspace_f, ntime;
  pk_matrix_->get_discretization_order(&nspace_m, &ntime);
  pk_fracture_->get_discretization_order(&nspace_f, &ntime);

  if (nspace_m == 2 && nspace_f == 2) {
    PK_MPCStrong<PK_BDF>::CalculateDiagnostics(tag);
  }
}


/* *******************************************************************
* Residual evaluation
******************************************************************* */
void TransportMatrixFractureImplicit_PK::FunctionalResidual(
     double t_old, double t_new,
     Teuchos::RCP<TreeVector> u_old,
     Teuchos::RCP<TreeVector> u_new,
     Teuchos::RCP<TreeVector> f)
{
  PK_MPCStrong<PK_BDF>::FunctionalResidual(t_old, t_new, u_old, u_new, f);

  // contribution from coupling terms
  auto u0 = u_new->SubVector(0)->Data();
  auto u1 = u_new->SubVector(1)->Data();

  auto f0 = f->SubVector(0)->Data();
  auto f1 = f->SubVector(1)->Data();

  op_coupling00_->Setup(fia_->get_values(), 1.0);
  op_coupling00_->global_operator()->Apply(*u0, *f0, 1.0);

  op_coupling01_->Setup(fia_->get_values2(), -1.0);
  op_coupling01_->global_operator()->Apply(*u1, *f0, 1.0);

  op_coupling10_->Setup(fia_->get_values(), -1.0);
  op_coupling10_->global_operator()->Apply(*u0, *f1, 1.0);

  op_coupling11_->Setup(fia_->get_values2(), 1.0);
  op_coupling11_->global_operator()->Apply(*u1, *f1, 1.0);

  if (flag_dispersion_) {
    const auto& kn = *S_->Get<CV_t>("fracture-normal_diffusion").ViewComponent("cell");
    fid_->SetValues(kn, 1.0);

    op_coupling00d_->Setup(fid_->get_values(), 1.0);
    op_coupling00d_->global_operator()->Apply(*u0, *f0, 1.0);

    op_coupling01d_->Setup(fid_->get_values(), -1.0);
    op_coupling01d_->global_operator()->Apply(*u1, *f0, 1.0);

    op_coupling10d_->Setup(fid_->get_values(), -1.0);
    op_coupling10d_->global_operator()->Apply(*u0, *f1, 1.0);

    op_coupling11d_->Setup(fid_->get_values(), 1.0);
    op_coupling11d_->global_operator()->Apply(*u1, *f1, 1.0);
  }
}


/* ******************************************************************
* Recompute preconditioner
****************************************************************** */
void TransportMatrixFractureImplicit_PK::UpdatePreconditioner(
    double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(tp, u, dtp);
}


/* ******************************************************************
* Apply block-diagonal preconditioner inv(B) * X.                                                 
****************************************************************** */
int TransportMatrixFractureImplicit_PK::ApplyPreconditioner(
    Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  return PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);
  // return op_pc_solver_->ApplyInverse(*X, *Y);
}

}  // namespace Amanzi

