#include "Reduced2p2c_PK.hh"
#include "EpetraExt_MultiVectorOut.h"
#include <EpetraExt_MatrixMatrix.h>

namespace Amanzi {
namespace Multiphase {

Reduced2p2c_PK::Reduced2p2c_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln):
                    soln_(soln), S_(S), pk_tree_(pk_tree), glist_(global_list),
                    passwd_("state")
{
  ti_list_ = Teuchos::rcp(new Teuchos::ParameterList(global_list->sublist("cycle driver").sublist("time integrator")));

  mpmc_list_ = Teuchos::sublist(global_list, "MPMC Specs", true);
  jacobian_type_ = mpmc_list_->get<std::string>("Jacobian type", "analytic");
  ncp_type_ = mpmc_list_->get<std::string>("NCP function", "min");
  mu_ = mpmc_list_->get<double>("smoothing parameter mu", 0.0);

  // We also need preconditioner and solver sublists
  pc_list_ = Teuchos::sublist(global_list, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(global_list, "solvers", true);

  // we need flow list
  cpr_enhanced_ = mpmc_list_->get<bool>("CPR enhancement", false);
  if (cpr_enhanced_) {
    cpr_list_ = Teuchos::sublist(mpmc_list_, "CPR parameters", true);
    std::vector<int> correction_blocks = cpr_list_->get<Teuchos::Array<int> >("correction blocks").toVector();
    pc_block_names_ = cpr_list_->get<Teuchos::Array<std::string> >("preconditioner").toVector();
    pc_block_names_.resize(correction_blocks.size());
    for (int i = 0; i < pc_block_names_.size(); i++)
      AMANZI_ASSERT(pc_list_->isSublist(pc_block_names_[i]));
  }

  pc_all_name_ = ti_list_->get<std::string>("preconditioner");
  linear_solver_name_ = ti_list_->get<std::string>("linear solver");
  AMANZI_ASSERT(pc_list_->isSublist(pc_all_name_));

  p1_tree_ = Teuchos::rcp(new TreeVector());
  s1_tree_ = Teuchos::rcp(new TreeVector());
  rhl_tree_ = Teuchos::rcp(new TreeVector());
  soln_->PushBack(p1_tree_);
  soln_->PushBack(s1_tree_);
  soln_->PushBack(rhl_tree_);
  comp_w_pk_ = Teuchos::rcp(new CompW_PK(pk_tree.sublist("Component 1"), mpmc_list_, S, soln_));
  comp_h_pk_ = Teuchos::rcp(new CompH_PK(pk_tree.sublist("Component 2"), mpmc_list_, S, soln_));
  gas_constraint_pk_ = Teuchos::rcp(new GasConstraint(*mpmc_list_, S));
  gas_constraint_pk_->setMu(mu_);
  // gas_constraint_pk_->SetNCPFunctionType(ncp_type_);

  num_mat_ = 0;
  ln_itrs_ = 0;
  nl_itrs_ = 0;
  ts_cnt_ = 0;

  // verbose object
  vo_ = new VerboseObject("Reduced2p2c::", *mpmc_list_);

  Amanzi::timer_manager.add("UpdatePreconditioner", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("ApplyPreconditioner", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("PreconditionerSetup", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("PreconditionerSolve", Amanzi::Timer::ACCUMULATE);

  accumulateSolveTime_ = 0.0;
  accumulateSetupTime_ = 0.0;
}


Reduced2p2c_PK::~Reduced2p2c_PK() {
  // Do nothing for now
}


void Reduced2p2c_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  rhs_ = Teuchos::rcp(new TreeVector());
  p1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("pressure_w",passwd_)));
  s1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("saturation_w",passwd_)));
  rhl_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("hydrogen density liquid",passwd_)));
  p1_tree_->SetData(p1_);
  s1_tree_->SetData(s1_);
  rhl_tree_->SetData(rhl_);
  comp_w_pk_->Initialize(S_.ptr());
  comp_h_pk_->Initialize(S_.ptr());
  gas_constraint_pk_->Initialize();

  auto cvs = Teuchos::rcp(new CompositeVectorSpace(comp_w_pk_->OpPrec1()->global_operator()->DomainMap()));
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> cvs_as_tvs = Teuchos::rcp(new TreeVectorSpace());
  cvs_as_tvs->SetData(cvs);
  tvs->PushBack(cvs_as_tvs);
  tvs->PushBack(cvs_as_tvs);
  tvs->PushBack(cvs_as_tvs);

  tree_op_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  tree_op_precond_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  if (cpr_enhanced_) 
    comb_tree_op_ = Teuchos::rcp(new Operators::CombinativeTreeOperator(tree_op_, cpr_list_, false));

  if( jacobian_type_ == "analytic") {
    tree_op_->SetOperatorBlock(0,0,comp_w_pk_->OpPrec1()->global_operator());
    tree_op_->SetOperatorBlock(0,1,comp_w_pk_->OpPrec2()->global_operator());
    tree_op_->SetOperatorBlock(0,2,comp_w_pk_->OpPrec3()->global_operator());
    tree_op_->SetOperatorBlock(1,0,comp_h_pk_->OpPrec1()->global_operator());
    tree_op_->SetOperatorBlock(1,1,comp_h_pk_->OpPrec2()->global_operator());
    tree_op_->SetOperatorBlock(1,2,comp_h_pk_->OpPrec3()->global_operator());

    tree_op_precond_->SetOperatorBlock(0,0,comp_w_pk_->OpPrec1()->global_operator());
    tree_op_precond_->SetOperatorBlock(0,1,comp_w_pk_->OpPrec2()->global_operator());
    tree_op_precond_->SetOperatorBlock(0,2,comp_w_pk_->OpPrec3()->global_operator());
    tree_op_precond_->SetOperatorBlock(1,0,comp_h_pk_->OpPrec1()->global_operator());
    tree_op_precond_->SetOperatorBlock(1,1,comp_h_pk_->OpPrec2()->global_operator());
    tree_op_precond_->SetOperatorBlock(1,2,comp_h_pk_->OpPrec3()->global_operator());
  } else if (jacobian_type_ == "numerical")
  {
    tree_op_->SetOperatorBlock(0,0,comp_w_pk_->Ops()[0]->global_operator());
    tree_op_->SetOperatorBlock(0,1,comp_w_pk_->Ops()[1]->global_operator());
    tree_op_->SetOperatorBlock(0,2,comp_w_pk_->Ops()[2]->global_operator());
    tree_op_->SetOperatorBlock(1,0,comp_h_pk_->Ops()[0]->global_operator());
    tree_op_->SetOperatorBlock(1,1,comp_h_pk_->Ops()[1]->global_operator());
    tree_op_->SetOperatorBlock(1,2,comp_h_pk_->Ops()[2]->global_operator());

    tree_op_precond_->SetOperatorBlock(0,0,comp_w_pk_->Ops()[0]->global_operator());
    tree_op_precond_->SetOperatorBlock(0,1,comp_w_pk_->Ops()[1]->global_operator());
    tree_op_precond_->SetOperatorBlock(0,2,comp_w_pk_->Ops()[2]->global_operator());
    tree_op_precond_->SetOperatorBlock(1,0,comp_h_pk_->Ops()[0]->global_operator());
    tree_op_precond_->SetOperatorBlock(1,1,comp_h_pk_->Ops()[1]->global_operator());
    tree_op_precond_->SetOperatorBlock(1,2,comp_h_pk_->Ops()[2]->global_operator());
  }
  tree_op_->SetOperatorBlock(2,0,gas_constraint_pk_->op_prec1()->global_operator());
  tree_op_->SetOperatorBlock(2,1,gas_constraint_pk_->op_prec2()->global_operator());
  tree_op_->SetOperatorBlock(2,2,gas_constraint_pk_->op_prec3()->global_operator());

  tree_op_precond_->SetOperatorBlock(2,0,gas_constraint_pk_->op_prec1_tmp()->global_operator());
  tree_op_precond_->SetOperatorBlock(2,1,gas_constraint_pk_->op_prec2_tmp()->global_operator());
  tree_op_precond_->SetOperatorBlock(2,2,gas_constraint_pk_->op_prec3_tmp()->global_operator());

  // Init time interval
  ProcessSublistTimeInterval(*ti_list_, ti_specs_generic_);
 
  ti_specs_generic_.T0  = ti_list_->get<double>("start interval time", 0);
  ti_specs_generic_.dT0 = ti_list_->get<double>("initial time step", 1);

  double T0 = ti_specs_generic_.T0;
  double dT0 = ti_specs_generic_.dT0;

  dT = dT0;
  dTnext = dT0;

  ti_specs_ = &ti_specs_generic_;

  error_control_ = ti_specs_->error_control_options;
  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
    ti_specs_->error_control_options = error_control_;
  }
  // set up new time integration or solver
  std::string ti_method_name(ti_specs_generic_.ti_method_name);

  if (ti_specs_generic_.ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList bdf1_list = ti_specs_generic_.ti_list_ptr_->sublist("BDF1");
    bdf1_dae = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  }

  // Initialize linear solvers
  if (cpr_enhanced_) 
    solver_comb_ = factory_comb.Create(linear_solver_name_, *linear_operator_list_, comb_tree_op_);
  else 
    solver_tree_ = factory_tree.Create(linear_solver_name_, *linear_operator_list_, tree_op_);

  // Initialize coarse indices array
  coarse_indices_array_ = new int[9];
  coarse_indices_array_[0] = 1;
  coarse_indices_array_[1] = 1;
  coarse_indices_array_[2] = 1;
  coarse_indices_array_[3] = 1;
  coarse_indices_array_[4] =-1;
  coarse_indices_array_[5] = 0;
  coarse_indices_array_[6] =-1;
  coarse_indices_array_[7] = 0;
  coarse_indices_array_[8] = 0;

  coarse_indices_array_two_level_ = new int[9];
  coarse_indices_array_two_level_[0] = 1;
  coarse_indices_array_two_level_[1] = 1;
  //coarse_indices_array_two_level_[2] = 1;
  coarse_indices_array_two_level_[3] = 1;
  coarse_indices_array_two_level_[4] = -1;
  // coarse_indices_array_two_level_[5] = 0;
  coarse_indices_array_two_level_[6] = -1;
  coarse_indices_array_two_level_[7] = 0;
  // coarse_indices_array_two_level_[8] = 0;

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << std::endl;
  std::string method_name = "Direct application of " + pc_all_name_;
  if (cpr_enhanced_) {
    if (pc_all_name_ != "Euclid")
      *vo_->os() << "Warning!!! Global Preconditioner is not Euclid. Solver may fail." << std::endl;
    std::stringstream ss;
    ss << pc_block_names_.size();   
    method_name = "CPR-AMG(" + ss.str() + ")";
  }

  *vo_->os() << "Jacobian type: " << jacobian_type_ << std::endl
    << "Preconditioning Method: " << method_name
    << std::endl << std::endl;
  *vo_->os() << "NCP function: " << ncp_type_ << std::endl;

  if (ncp_type_ == "fischer-burmeister")  
    *vo_->os() << "smoothing parameter mu: " << mu_ << std::endl;

  *vo_->os() << vo_->color("green") << "Initializaton complete" << vo_->reset() << std::endl;
}


bool Reduced2p2c_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dT = t_new - t_old;
  dT_actual = dT;
  double state_time = S_->time();
  double time = 0.0;
  if (time >= 0.0) T_physics = time;
  time = T_physics;
  // int cycle = S_->cycle();
  // if (cycle == 100) comp_h_pk_->ComputeBCs(cycle);

  double time_stop_injection = 1.5768e+13;
  if (state_time - time_stop_injection > 1e-10) {
    comp_h_pk_->ComputeBCs(true);
  }

  Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(soln_->Map()));
  udot->PutScalar(0.0);

  /*
  // predict water mass change during time step
  if (ti_specs_->num_itrs == 0) {  // initialization
    UpdatePreconditioner(time, soln_, dT);
    ti_specs_->num_itrs++;
  }
  */

  bdf1_dae->SetInitialState(time, soln_, udot);

  bool fail = false;
  if (ti_specs_->ti_method == FLOW_TIME_INTEGRATION_BDF1){
    fail = bdf1_dae->TimeStep(dT, dTnext, soln_);
    if (fail) {
      dT = dTnext;
      nl_itrs_ = 0;
      return fail;
    }
  }

  // tell the caller what time step we actually took
  dT_actual = dT;
  
  nl_itrs_ = 0;
  ti_specs_->num_itrs++;
  dT = dTnext;

  return fail;
}


void Reduced2p2c_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  double mu_k;
  if (nl_itrs_ == 0) ts_cnt_++;
  if (nl_itrs_ <= 4) {
    mu_k = pow(0.1, nl_itrs_) * mu_;
    gas_constraint_pk_->setMu(mu_k);
  } else {
    mu_k = 0.0;
    gas_constraint_pk_->setMu(mu_k);
  }
  nl_itrs_++;
  Amanzi::timer_manager.start("UpdatePreconditioner");

  if (jacobian_type_ == "analytic") {
    comp_w_pk_->UpdatePreconditioner(t, up, h);
    comp_h_pk_->UpdatePreconditioner(t, up, h);
  } else {
    // comp_w_pk_->NumericalJacobian(0.0, h, up, 1e-12);
    // comp_h_pk_->NumericalJacobian(0.0, h, up, 1e-12);
  }
  gas_constraint_pk_->UpdatePreconditioner(t, up, h);

  if (pc_list_->sublist(pc_all_name_).get<std::string>("preconditioner type") == "systg") {
    Teuchos::ParameterList& systg_list =  pc_list_->sublist(pc_all_name_).sublist("systg parameters");
    std::cout << "mu_k = " << mu_k << std::endl;
    if (std::abs(mu_k) < 1e-15) {
      systg_list.set("max coarse levels", 3);
      systg_list.set("inactive gas indices", gas_constraint_pk_->getInactiveGasIndices());
      systg_list.set("number of inactive cells", gas_constraint_pk_->getNumInactiveCells());
      systg_list.set("coarse indices array", coarse_indices_array_);
    } else {
      systg_list.set("max coarse levels", 2);
      systg_list.set("inactive gas indices", gas_constraint_pk_->getInactiveGasIndices());
      systg_list.set("number of inactive cells", 0);
      systg_list.set("coarse indices array", coarse_indices_array_two_level_);
    }
  }

  if (cpr_enhanced_) {
    comb_tree_op_->SymbolicAssembleMatrix();
    comb_tree_op_->AssembleMatrix();
    //comb_tree_op_->InitPreconditionerGlobal(pc_all_name_, *pc_list_);
    //comb_tree_op_->InitPreconditionerBlock(pc_block_names_, *pc_list_);   
    //ComputeAPinv(); 
  } else {
    tree_op_->SymbolicAssembleMatrix();
    tree_op_->AssembleMatrix();
    //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());
    /*
    if (gas_constraint_pk_->getNumCells() < 2000) {
      Teuchos::RCP<Epetra_Vector> diag_values = Teuchos::rcp(new Epetra_Vector(tree_op_->A()->RowMap()));
      tree_op_->A()->ExtractDiagonalCopy(*diag_values);
      Teuchos::RCP<Epetra_Vector> shift_values = Teuchos::rcp(new Epetra_Vector(*diag_values));
      shift_values->PutScalar(0.0);
      for (int i = 0; i < shift_values->MyLength(); i++) {
        if ((i % 3) == 0) {
          (*shift_values)[i] = 1.0e-3;
        }
      }
      diag_values->Update(1.0, *shift_values, 1.0);
      int err = tree_op_->A()->ReplaceDiagonalValues(*diag_values);
      ASSERT(err == 0);
    }
    */
    //tree_op_->InitPreconditioner(pc_all_name_, *pc_list_);
  }

  // tree_op_->SymbolicAssembleMatrix();
  // tree_op_->AssembleMatrix();
  Amanzi::timer_manager.stop("UpdatePreconditioner");
}


void Reduced2p2c_PK::FunctionalResidual(double t_old, double t_new, 
                                        Teuchos::RCP<TreeVector> u_old,
                                        Teuchos::RCP<TreeVector> u_new,
                                        Teuchos::RCP<TreeVector> f) 
{
  comp_w_pk_->FunctionalResidual(t_old, t_new, u_old, u_new, f->SubVector(0));
  comp_h_pk_->FunctionalResidual(t_old, t_new, u_old, u_new, f->SubVector(1));
  gas_constraint_pk_->FunctionalResidual(t_old, t_new, u_old, u_new, f->SubVector(2));
  rhs_ = f;
}


int Reduced2p2c_PK::ApplyJacobian(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu)
{
  tree_op_->Apply(*u, *pu);
  return 0;
}


int Reduced2p2c_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu)
{
  Amanzi::timer_manager.start("ApplyPreconditioner");

  /*if (ln_itrs_ == 0)*/ pu->PutScalar(0.0);
  TreeVector pu_temp = *pu;
  num_mat_++;
  std::stringstream ss;
  ss << num_mat_;
  std::string file_name = "matrixA_" + ss.str() + ".txt";
  int ierr = 0;
  if (cpr_enhanced_) {
    comb_tree_op_->InitPreconditionerGlobal(pc_all_name_, *pc_list_);
    comb_tree_op_->InitPreconditionerBlock(pc_block_names_, *pc_list_);  
    ierr = solver_comb_->ApplyInverse(*u, *pu);
    if (ierr > 0) ln_itrs_ += solver_comb_->num_itrs();
    //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *comb_tree_op_->Op()->A());
  }
  else {
    //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());
    Amanzi::timer_manager.start("PreconditionerSetup");
    tree_op_->InitPreconditioner(pc_all_name_, *pc_list_);
    Amanzi::timer_manager.stop("PreconditionerSetup");
    Amanzi::timer_manager.start("PreconditionerSolve");
    ierr = solver_tree_->ApplyInverse(*u, *pu);
    Amanzi::timer_manager.stop("PreconditionerSolve");
    if (ierr > 0) 
      ln_itrs_ += solver_tree_->num_itrs();
    else {
      EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());
      std::ofstream file_out;
      file_name = "rhs_" + ss.str() + ".txt";
      file_out.open(file_name.c_str());
      int ncells = (*u->SubVector(0)->Data()->ViewComponent("cell")).MyLength();
      const Epetra_MultiVector& rhs_p = *u->SubVector(0)->Data()->ViewComponent("cell");
      const Epetra_MultiVector& rhs_s = *u->SubVector(1)->Data()->ViewComponent("cell");
      const Epetra_MultiVector& rhs_c = *u->SubVector(2)->Data()->ViewComponent("cell");
      for (int i = 0; i < ncells; i++) {
        file_out << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << rhs_p[0][i] << std::endl 
                 << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << rhs_s[0][i] << std::endl 
                 << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << rhs_c[0][i] << std::endl;
      }
      file_out.close();

      file_name = "inactiveIdx_" + ss.str() + ".txt";
      file_out.open(file_name.c_str());
      file_out << gas_constraint_pk_->getNumInactiveCells() << std::endl;
      for (int i=0; i < gas_constraint_pk_->getNumInactiveCells(); i++) {
        file_out << *(gas_constraint_pk_->getInactiveGasIndices() + i) << std::endl;
      }
      file_out.close();
    }
  }

  if (ierr < 0) {
    Teuchos::ParameterList& systg_list =  pc_list_->sublist(pc_all_name_).sublist("systg parameters");
    systg_list.set("coarse solver type", 1);
    systg_list.set("coarse solver iter", 5);
    *vo_->os() << vo_->color("green") << "Use GMRES for coarse solver\n";
    ierr = 0;
    tree_op_->InitPreconditioner(pc_all_name_, *pc_list_);
    ierr = solver_tree_->ApplyInverse(*u, pu_temp);
    *pu = pu_temp;
    if (ierr > 0) ln_itrs_ += solver_tree_->num_itrs();
  }

  //AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> factory;
  //Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::TreeOperator, TreeVector, TreeVectorSpace> >
  //    solver = factory.Create("GMRES", *linear_operator_list_, tree_op_);
  //solver->ApplyInverse(*u, *pu);
  //ln_itrs_ += solver->num_itrs();

  //tree_op_->ApplyInverse(*u, *pu);

  /*
  Teuchos::RCP<TreeVector> out = Teuchos::rcp(new TreeVector(*u));
  out->PutScalar(0.0);
  tree_op_->Apply(*pu, *out);
  */

  Amanzi::timer_manager.stop("ApplyPreconditioner");
  return ierr;
}


double Reduced2p2c_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) 
{
  // Monitor l2 norm of residual
  double du_l2 = 0.0;
  double resnorm_p, resnorm_s, resnorm_r;
  du->SubVector(0)->Data()->Norm2(&resnorm_p);
  du->SubVector(1)->Data()->Norm2(&resnorm_s);
  du->SubVector(2)->Data()->Norm2(&resnorm_r);
  printf("resnorm_p = %4.6e, resnorm_s = %4.6e, resnorm_r = %4.6e \n", resnorm_p, resnorm_s, resnorm_r);
  du->Norm2(&du_l2);
  return du_l2;
}


void Reduced2p2c_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  comp_w_pk_->CommitStep(t_old, t_new, S);
  comp_h_pk_->CommitStep(t_old, t_new, S);
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Reduced2p2c_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
           Teuchos::RCP<const TreeVector> u,
           Teuchos::RCP<TreeVector> du)
{
  // du->Scale(h);
  /*
  Teuchos::RCP<CompositeVector> s_next = Teuchos::rcp(new CompositeVector(*u->SubVector(1)->Data()));
  s_next->Update(-1.0, *du->SubVector(1)->Data(), 1.0);
  ClipSaturation(s_next, 1e-12);
  du->SubVector(1)->Data()->Update(1.0, *u->SubVector(1)->Data(), -1.0, *s_next, 0.0);
  */
  Teuchos::RCP<CompositeVector> rho_next = Teuchos::rcp(new CompositeVector(*u->SubVector(2)->Data()));
  rho_next->Update(-1.0, *du->SubVector(2)->Data(), 1.0);
  ClipConcentration(rho_next);
  du->SubVector(2)->Data()->Update(1.0, *u->SubVector(2)->Data(), -1.0, *rho_next, 0.0);

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}


void Reduced2p2c_PK::ClipSaturation(Teuchos::RCP<CompositeVector> s, double tol)
{
  double max_s, min_s;
  max_s = 1.0 - tol;
  min_s = 0.0 + tol;
  Epetra_MultiVector& s_c = *s->ViewComponent("cell");
  for (int c = 0; c < s_c.MyLength(); c++)
  {
    s_c[0][c] = std::min(std::max(min_s, s_c[0][c]), max_s);
  }
}

void Reduced2p2c_PK::ClipConcentration(Teuchos::RCP<CompositeVector> rho)
{
  Epetra_MultiVector& rho_c = *rho->ViewComponent("cell");
  for (int c = 0; c < rho_c.MyLength(); c++)
  {
    rho_c[0][c] = std::max(0.0, rho_c[0][c]);
  }
}


/* ******************************************************************
* Process generic time  interval sublist.
**************************************************************** */
void Reduced2p2c_PK::ProcessSublistTimeInterval(
  Teuchos::ParameterList& ti_list,  TI_Specs& ti_specs){

  ti_specs.ti_list_ptr_ = &ti_list;
  std::string ti_method_name = ti_list.get<std::string>("time integration method", "none");

  ProcessStringTimeIntegration(ti_method_name, &ti_specs.ti_method);
  ProcessSublistTimeIntegration(ti_list, ti_method_name, ti_specs);
  ti_specs.ti_method_name = ti_list.name();

  ti_specs.preconditioner_name = FindStringPreconditioner(ti_list);
  ti_specs.solver_name = FindStringLinearSolver(ti_list);

  ProcessStringErrorOptions(ti_list, &ti_specs.error_control_options);
}


/* ******************************************************************
* Process time integration sublist.
**************************************************************** */
void Reduced2p2c_PK::ProcessSublistTimeIntegration(
    Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs)
{
  Errors::Message msg;

  if (list.isSublist(name)) {
    Teuchos::ParameterList bdf1_list = list.sublist(name);

    // new way to define parameters overrides the above values.
    if (bdf1_list.isParameter("timestep controller type")) {
      std::string dT_method_name = bdf1_list.get<std::string>("timestep controller type");

      ti_specs.dT_method = 0;
      Teuchos::ParameterList dtlist;
      if (dT_method_name == "standard") {
        dtlist = bdf1_list.sublist("timestep controller standard parameters");
        ti_specs.dTfactor = dtlist.get<double>("time step increase factor");
      } else if (dT_method_name == "fixed") {
        dtlist = bdf1_list.sublist("timestep controller fixed parameters");
        ti_specs.dTfactor = dtlist.get<double>("time step increase factor");
      } else if (dT_method_name == "adaptive") {
        dtlist = bdf1_list.sublist("timestep controller adaptive parameters");
        ti_specs.dT_method = FLOW_DT_ADAPTIVE;
      }
      ti_specs.dTmax = dtlist.get<double>("max time step", FLOW_MAXIMUM_DT);
    }

    if (list.isSublist("initialization")) {
      Teuchos::ParameterList& ini_list = list.sublist("initialization");
      std::string name = ini_list.get<std::string>("method", "none");
      ti_specs.initialize_with_darcy = (name == "saturated solver");
      ti_specs.clip_saturation = ini_list.get<double>("clipping saturation value", -1.0);
      ti_specs.clip_pressure = ini_list.get<double>("clipping pressure value", -1e+10);

      ti_specs.solver_name_ini = FindStringLinearSolver(ini_list);
    }

    if (list.isSublist("pressure-lambda constraints")) {
      Teuchos::ParameterList& pl_list = list.sublist("pressure-lambda constraints");
      ti_specs.pressure_lambda_constraints = true;
      ti_specs.inflow_krel_correction = pl_list.get<bool>("inflow krel correction", "false");

      ti_specs.solver_name_constraint = FindStringLinearSolver(pl_list);
    } else {
      ti_specs.pressure_lambda_constraints = false;
    }

    // Picard sublist
    ti_specs.max_itrs = bdf1_list.get<int>("maximum number of iterations", FLOW_TI_MAX_ITERATIONS);
    ti_specs.residual_tol = bdf1_list.get<double>("convergence tolerance", FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE);

  } else if (name != "none") {
    msg << "\n MPMC PK: specified time integration sublist does not exist.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the time integration method.
**************************************************************** */
void Reduced2p2c_PK::ProcessStringTimeIntegration(const std::string name, int* method)
{
  Errors::Message msg;

  if (name == "Picard") {
    *method = FLOW_TIME_INTEGRATION_PICARD;
  } else if (name == "backward Euler") {
    *method = FLOW_TIME_INTEGRATION_BACKWARD_EULER;
  } else if (name == "BDF1") {
    *method = FLOW_TIME_INTEGRATION_BDF1;
  } else if (name == "BDF2") {
    *method = FLOW_TIME_INTEGRATION_BDF2;
  } else {
    msg << "MPMC PK: time integration method \"" << name.c_str() << "\" is not known.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
std::string Reduced2p2c_PK::FindStringPreconditioner(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("preconditioner")) {
    name = list.get<std::string>("preconditioner");
  } else {
    msg << "MPMC PK: parameter <preconditioner> is missing either in TI or LS list.";
    Exceptions::amanzi_throw(msg);
  }

  if (! pc_list_->isSublist(name)) {
    msg << "MPMC PK: preconditioner \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Find string for the linear solver
**************************************************************** */
std::string Reduced2p2c_PK::FindStringLinearSolver(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("linear solver")) {
    name = list.get<std::string>("linear solver");
  } else {
    msg << "MPMC PK: time integrator does not define <linear solver>.";
    Exceptions::amanzi_throw(msg);
  }

  if (! linear_operator_list_->isSublist(name)) {
    msg << "MPMC PK: linear solver \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Process string for error control options
**************************************************************** */
void Reduced2p2c_PK::ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control)
{
  *control = 0;
  if (list.isParameter("error control options")){
    std::vector<std::string> options;
    options = list.get<Teuchos::Array<std::string> >("error control options").toVector();

    for (int i=0; i < options.size(); i++) {
      if (options[i] == "pressure") {
        *control += FLOW_TI_ERROR_CONTROL_PRESSURE;
      } else if (options[i] == "saturation") {
        *control += FLOW_TI_ERROR_CONTROL_SATURATION;
      } else if (options[i] == "residual") {
        *control += FLOW_TI_ERROR_CONTROL_RESIDUAL;
      } else {
        Errors::Message msg;
        msg << "MPMC PK: unknown error control option has been specified.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
}


}
}
