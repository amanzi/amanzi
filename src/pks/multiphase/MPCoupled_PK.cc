#include "MPCoupled_PK.hh"
#include "EpetraExt_MultiVectorOut.h"
#include <EpetraExt_MatrixMatrix.h>

namespace Amanzi {
namespace Multiphase {

MPCoupled_PK::MPCoupled_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln):
                    soln_(soln), S_(S), pk_tree_(pk_tree), glist_(global_list),
                    passwd_("state"), include_capillary_(false), cpr_enhanced_(true),
                    include_gravity_(false)
{
  // We also need preconditioner and solver sublists
  pc_list_ = Teuchos::sublist(global_list, "Preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(global_list, "Solvers", true);

  // we need flow list
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(global_list, "Flow", true);
  jacobian_type_ = flow_list->get<std::string>("Jacobian type", "analytic");
  include_capillary_ = flow_list->get<bool>("Include capillary pressure", false);
  use_total_mobility_ = flow_list->get<bool>("use total mobility", false);
  include_gravity_ = 
    flow_list->sublist("Phase1 Specs").sublist("operators").sublist("diffusion operator").sublist("matrix").get<bool>("gravity",false) &&
    flow_list->sublist("Phase2 Specs").sublist("operators").sublist("diffusion operator").sublist("matrix").get<bool>("gravity",false);
  cpr_enhanced_ = flow_list->get<bool>("CPR enhancement", true);
  block_factorization_ = flow_list->get<bool>("use block factorization", false);
  cpr_list_ = Teuchos::sublist(flow_list, "CPR parameters", true);
  if (cpr_enhanced_) {
    std::vector<int> correction_blocks = cpr_list_->get<Teuchos::Array<int> >("correction blocks").toVector();
    pc_block_names_ = cpr_list_->get<Teuchos::Array<std::string> >("preconditioner").toVector();
    pc_block_names_.resize(correction_blocks.size());
    for (int i = 0; i < pc_block_names_.size(); i++)
      ASSERT(pc_list_->isSublist(pc_block_names_[i]));
  }

  ti_list_ = Teuchos::sublist(flow_list, "time integrator", true);
  pc_all_name_ = ti_list_->get<std::string>("preconditioner");
  pc_ass_name_ = flow_list->get<std::string>("A_ss preconditioner", pc_all_name_);
  linear_solver_name_ = ti_list_->get<std::string>("linear solver");
  ASSERT(pc_list_->isSublist(pc_all_name_));

  Teuchos::RCP<Teuchos::ParameterList> misc_list = Teuchos::sublist(global_list, "Misc", false);
  block01_scaling_ = misc_list->get<double>("scaling block01", 1.0);
  block10_scaling_ = misc_list->get<double>("scaling block10", 1.0);

  // create solution tree vector
  p1_tree_ = Teuchos::rcp(new TreeVector());
  s2_tree_ = Teuchos::rcp(new TreeVector());
  soln_->PushBack(p1_tree_);
  soln_->PushBack(s2_tree_);
  ln_itrs_ = 0;
  nl_itrs_ = 0;

  // create the pks for phases
  phase2_pk_ = Teuchos::rcp(new Phase2_PK(pk_tree, global_list, S, soln_));
  phase2_pk_->IncludeCapillary(include_capillary_);
  if (use_total_mobility_) {
    tot_phase_pk_ = Teuchos::rcp(new TotalPhase_PK(pk_tree, global_list, S, soln_));
    tot_phase_pk_->IncludeCapillary(include_capillary_);
  } else {
    phase1_pk_ = Teuchos::rcp(new Phase1_PK(pk_tree, global_list, S, soln_));
  }
  //std::cout << "Done MPCoupled_PK Constructor\n";

  // temporary counter for output matrices
  num_mat_ = 0;

  // verbose object
  vo_ = new VerboseObject("MPCoupled_PK::", *flow_list); 

  // Timer for profiling
  Amanzi::timer_manager.add("Functional", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("UpdatePreconditioner", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("ApplyPreconditioner", Amanzi::Timer::ACCUMULATE);
}

MPCoupled_PK::~MPCoupled_PK() {
  // Do nothing for now
}

void MPCoupled_PK::Initialize() {
  p1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("phase1_pressure")));
  s2_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("phase2_saturation")));

  p1_tree_->SetData(p1_);
  s2_tree_->SetData(s2_);
  if (use_total_mobility_) {
    tot_phase_pk_->Initialize();
  } else {
    phase1_pk_->Initialize();
  }

  phase2_pk_->Initialize();

  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> cvs_as_tv = Teuchos::rcp(new TreeVectorSpace(cvs));
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);

  // create the tree operator for the blocks obtained from phase_pks
  tree_op_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  if (cpr_enhanced_ || block_factorization_) 
    comb_tree_op_ = Teuchos::rcp(new Operators::CombinativeTreeOperator(tree_op_, cpr_list_, block_factorization_));

  if (jacobian_type_ == "analytic")
  {
    if (use_total_mobility_) {
      tree_op_->SetOperatorBlock(0,0,tot_phase_pk_->op_prec1()->global_operator());
      tree_op_->SetOperatorBlock(0,1,tot_phase_pk_->op_prec2()->global_operator());
    } else {
      tree_op_->SetOperatorBlock(0,0,phase1_pk_->op_prec1()->global_operator());
      tree_op_->SetOperatorBlock(0,1,phase1_pk_->op_prec2()->global_operator());     
    }
    tree_op_->SetOperatorBlock(1,0,phase2_pk_->op_prec1()->global_operator());
    tree_op_->SetOperatorBlock(1,1,phase2_pk_->op_prec2()->global_operator());
  } else if (jacobian_type_ == "numerical")
  {
    if (use_total_mobility_) {
      tree_op_->SetOperatorBlock(0,0,tot_phase_pk_->op_prec1()->global_operator());
      tree_op_->SetOperatorBlock(0,1,tot_phase_pk_->op_prec2()->global_operator());
    } else {
      tree_op_->SetOperatorBlock(0,0,phase1_pk_->op_prec1()->global_operator());
      tree_op_->SetOperatorBlock(0,1,phase1_pk_->op_prec2()->global_operator());     
    }
    tree_op_->SetOperatorBlock(1,0,phase2_pk_->Ops()[0]->global_operator());
    tree_op_->SetOperatorBlock(1,1,phase2_pk_->Ops()[1]->global_operator());
  }

  dt_ = ti_list_->get<double>("initial time step");
  dTnext = dt_;

  // set up new time integration or solver
  std::string ti_method_name = ti_list_->get<std::string>("time integration method");
  ASSERT(ti_method_name == "BDF1");

  Teuchos::RCP<Teuchos::ParameterList> bdf1_list = Teuchos::sublist(ti_list_, "BDF1", true);
  //std::cout<<bdf1_list<<"\n";
  //if (! bdf1_list.isSublist("VerboseObject"))
  //    bdf1_list.sublist("VerboseObject") = mp_list_.sublist("VerboseObject");
  bdf1_dae = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, *bdf1_list, soln_));

  if (cpr_enhanced_ || block_factorization_) 
    solver_comb_ = factory_comb.Create(linear_solver_name_, *linear_operator_list_, comb_tree_op_);
  else 
    solver_tree_ = factory_tree.Create(linear_solver_name_, *linear_operator_list_, tree_op_);

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << std::endl;
  std::string method_name = "Direct application of " + pc_all_name_;
  if (block_factorization_) method_name = "Block Factorization with " + pc_all_name_;
  if (cpr_enhanced_) {
    if (pc_all_name_ != "Euclid")
      *vo_->os() << "Warning!!! Global Preconditioner is not Euclid. Solver may fail." << std::endl;
    std::stringstream ss;
    ss << pc_block_names_.size();   
    method_name = "CPR-AMG(" + ss.str() + ")";
  }

  if (cpr_enhanced_ && block_factorization_) {
    *vo_->os() << "Warning!!! Cannot use CPR-AMG and Block Factorization at the same time." << std::endl
      << "Use Block Factorization as default." << std::endl;
    cpr_enhanced_ = false;
    method_name = "Block Factorization with " + pc_all_name_;
  }
    *vo_->os() << vo_->color("green") << "Initializaton Complete!"
      << vo_->reset() << std::endl;
    *vo_->os() << "Jacobian type: " << jacobian_type_ << std::endl
      << "Preconditioning Method: " << method_name
      << std::endl << std::endl;

  // Initialize coarse indices array
  coarse_indices_array_ = new int[2];
  coarse_indices_array_[0] = 1;
  coarse_indices_array_[1] = 0;
}


void MPCoupled_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  //std::cout << "UpdatePreconditioner: solution \n";
  //up->Print(std::cout);

  if (pc_list_->sublist(pc_all_name_).get<std::string>("preconditioner type") == "systg") {
    Teuchos::ParameterList& systg_list =  pc_list_->sublist(pc_all_name_).sublist("systg parameters");
    systg_list.set("coarse indices array", coarse_indices_array_);
    //std::cout << systg_list << "\n";
  }

  num_mat_++;
  Amanzi::timer_manager.start("UpdatePreconditioner");
  // we get the blocks from each pk
  if (jacobian_type_ == "analytic") {
    if (use_total_mobility_) {
      tot_phase_pk_->UpdatePreconditioner(t, up, h);
    } else {
      phase1_pk_->UpdatePreconditioner(t, up, h);
    }
    phase2_pk_->UpdatePreconditioner(t, up, h);
    //phase1_pk_->op_prec2()->global_operator()->Rescale(block01_scaling_);
    //phase2_pk_->op_prec1()->global_operator()->Rescale(block10_scaling_);
  } else if (jacobian_type_ == "numerical") {
    if (use_total_mobility_) {
      tot_phase_pk_->UpdatePreconditioner(t, up, h);
    } else {
      phase1_pk_->UpdatePreconditioner(t, up, h);
    }
    phase2_pk_->NumericalJacobian(0.0, h, up, 1e-12);
  }

  if (!block_factorization_) {
    if (cpr_enhanced_) {
      comb_tree_op_->SymbolicAssembleMatrix();
      comb_tree_op_->AssembleMatrix();
      //std::cout << "full matrix: " << *comb_tree_op_->Op()->A() << "\n";
      std::stringstream ss;
      ss << num_mat_;
      std::string file_name = "matrixA_" + ss.str() + ".txt";
      //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *comb_tree_op_->Op()->A());
      comb_tree_op_->InitPreconditionerGlobal(pc_all_name_, *pc_list_);
      comb_tree_op_->InitPreconditionerBlock(pc_block_names_, *pc_list_);   
      //ComputeAPinv(); 
    } else {
      tree_op_->SymbolicAssembleMatrix();
      tree_op_->AssembleMatrix();
      //std::cout << "full matrix: " << *tree_op_->A() << "\n";
      std::stringstream ss;
      ss << num_mat_;
      std::string file_name = "matrixA_" + ss.str() + ".txt";
      //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());
      tree_op_->InitPreconditioner(pc_all_name_, *pc_list_);
      //ComputeAPinv();
    }
  }

  /***************************************************
  *        Testing Schur Complement method           *
  ****************************************************/
  else {
    //std::cout << "UpdatePreconditioner\n";
    std::stringstream ss;
    ss << num_mat_;
    std::string file_name = "matrixA_" + ss.str() + ".txt";
    //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());

    //TreeVector unit_tv = *up;
    //unit_tv.PutScalar(1.0);
    //TreeVector tmp_res = unit_tv;
    //tmp_res.PutScalar(0.0);
    //tree_op_->Apply(unit_tv, tmp_res);
    //std::cout << "tmp_res before modification: \n";
    //tmp_res.Print(std::cout);

    // Assemble the all the blocks
    tree_op_->GetBlock(0,0)->SymbolicAssembleMatrix();
    tree_op_->GetBlock(0,0)->AssembleMatrix();
    tree_op_->GetBlock(1,1)->SymbolicAssembleMatrix();
    tree_op_->GetBlock(1,1)->AssembleMatrix();
    tree_op_->GetBlock(0,1)->SymbolicAssembleMatrix();
    tree_op_->GetBlock(0,1)->AssembleMatrix();
    tree_op_->GetBlock(1,0)->SymbolicAssembleMatrix();
    tree_op_->GetBlock(1,0)->AssembleMatrix();

    // extract the diagonal of A_ss and invert it
    Teuchos::RCP<Epetra_Vector> A_ss_diag = 
      Teuchos::rcp(new Epetra_Vector(tree_op_->GetBlock(1,1)->A()->RowMap()));
    tree_op_->GetBlock(1,1)->A()->ExtractDiagonalCopy(*A_ss_diag);
    //std::cout << *A_ss_diag << "\n";
    A_ss_diag->Reciprocal(*A_ss_diag);
    //std::cout << "inv diag: " << *A_ss_diag << "\n";

    Teuchos::RCP<Epetra_Vector> A_ps_diag = 
      Teuchos::rcp(new Epetra_Vector(tree_op_->GetBlock(0,1)->A()->RowMap()));
    tree_op_->GetBlock(0,1)->A()->ExtractDiagonalCopy(*A_ss_diag); 

    Teuchos::RCP<Epetra_CrsMatrix> A_sp = 
      Teuchos::rcp(new Epetra_CrsMatrix(*tree_op_->GetBlock(1,0)->A()));

    // create a Schur complement matrix from scratch
    // schur_mat should not be filled until all the operations are performed
    Teuchos::RCP<Epetra_CrsMatrix> schur_mat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, tree_op_->GetBlock(0,1)->A()->DomainMap(),1));
    schur_mat->PutScalar(0.0);
    // compute A_ps * inv(diag(A_ss)) * A_sp
    //tree_op_->GetBlock(1,0)->A()->LeftScale(*A_ss_diag);
    A_sp->LeftScale(*A_ss_diag);
    // use full A_ps
    //int err = EpetraExt::MatrixMatrix::Multiply(*tree_op_->GetBlock(0,1)->A(), false,
    //                                            *tree_op_->GetBlock(1,0)->A(), false, *schur_mat, false);
    //ASSERT(err == 0);

    // use diag(A_ps)
    A_sp->LeftScale(*A_ps_diag);
    int err = EpetraExt::MatrixMatrix::Add(*tree_op_->GetBlock(1,0)->A(), false, 1.0, *schur_mat, 0.0);
    ASSERT(err == 0);

    //std::cout << "B * inv(D) * C: \n" << *schur_mat << "\n";
    // S = A_pp - A_ps * inv(diag(A_ss)) * A_sp
    err = EpetraExt::MatrixMatrix::Add(*tree_op_->GetBlock(0,0)->A(), false, 1.0, *schur_mat, -1.0);
    ASSERT(err == 0);
    // now fill the matrix and modify the operator
    schur_mat->FillComplete();
    Teuchos::RCP<Epetra_CrsMatrix> A_pp = tree_op_->GetBlock(0,0)->A();
    tree_op_->GetBlock(0,0)->A() = schur_mat;

    //file_name = "schur_" + ss.str() + ".txt";
    //EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->GetBlock(0,0)->A());
    //tree_op_->Apply(unit_tv, tmp_res);
    //std::cout << "tmp_res after modification: \n";
    //tmp_res.Print(std::cout);

    tree_op_->GetBlock(0,0)->InitPreconditioner(pc_all_name_, *pc_list_);
    tree_op_->GetBlock(0,0)->A() = A_pp;
    //std::cout << "AFF preconditioner: " << pc_ass_name_ << "\n";
    tree_op_->GetBlock(1,1)->InitPreconditioner(pc_ass_name_, *pc_list_);    
    comb_tree_op_->SymbolicAssembleMatrix();
    comb_tree_op_->AssembleMatrix();
  }

  Amanzi::timer_manager.stop("UpdatePreconditioner");
  //MPI_Comm mpi_comm(MPI_COMM_WORLD);
  //Amanzi::timer_manager.parSync(mpi_comm);
  //Amanzi::timer_manager.print();
}


void MPCoupled_PK::Functional(double t_old, double t_new, 
                         Teuchos::RCP<TreeVector> u_old,
                         Teuchos::RCP<TreeVector> u_new,
                         Teuchos::RCP<TreeVector> f) 
{
  Amanzi::timer_manager.start("Functional");
  //std::cout << "u_new: \n";
  //u_new->Print(std::cout);
  if (use_total_mobility_) {
    tot_phase_pk_->Functional(t_old, t_new, u_old, u_new, f->SubVector(0));
  } else {
    phase1_pk_->Functional(t_old, t_new, u_old, u_new, f->SubVector(0));
  }
  phase2_pk_->Functional(t_old, t_new, u_old, u_new, f->SubVector(1));
  //std::cout << "functional phase1: " << *f->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  //std::cout << "functional phase2: " << *f->SubVector(1)->Data()->ViewComponent("cell") << "\n";

  /*
  ofstream rhs;
  rhs.open("rhs.txt");
  int ncells = (*f->SubVector(0)->Data()->ViewComponent("cell")).MyLength();
  const Epetra_MultiVector& rhs_p = *f->SubVector(0)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& rhs_s = *f->SubVector(1)->Data()->ViewComponent("cell");
  for (int i = 0; i < ncells; i++) {
    rhs << rhs_p[0][i] << endl << rhs_s[0][i] << endl;
  }
  rhs.close();
  */
  Amanzi::timer_manager.stop("Functional");
}


bool MPCoupled_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  dt_ = t_new - t_old;
  double time = S_->time();
  if (time >= 0.0) T_physics = time;
  time = T_physics;

  if (include_gravity_) {
    std::cout << "Make 1 iteration\n";
    // make one iteration for pressure
    //Teuchos::RCP<TreeVector> tmp_pres = Teuchos::rcp(new TreeVector(*soln_->SubVector(0)));
    Teuchos::RCP<TreeVector> tmp_pres_res = Teuchos::rcp(new TreeVector(*soln_->SubVector(0)));
    UpdatePreconditioner(0.0, soln_, 0.0);
    phase1_pk_->Functional(t_old, t_new, soln_, soln_, tmp_pres_res);
    if (nl_itrs_ == 0) {
      tmp_pres_res->Scale(-1.0);
      if (!cpr_enhanced_) {
        tree_op_->GetBlock(0,0)->SymbolicAssembleMatrix();
        tree_op_->GetBlock(0,0)->AssembleMatrix();
        tree_op_->GetBlock(0,0)->InitPreconditioner(pc_all_name_, *pc_list_);
      }
      tree_op_->GetBlock(0,0)->ApplyInverse(*tmp_pres_res->Data(), *soln_->SubVector(0)->Data());
    }
    //std::cout << "pressure after 1 iter: \n";
    //soln_->SubVector(0)->Data()->Print(std::cout);
  }

  Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(soln_->Map()));
  udot->PutScalar(0.0);

  bdf1_dae->SetInitialState(time, soln_, udot);
  bool fail = false;
  fail = bdf1_dae->TimeStep(dt_, dTnext, soln_);
  if (fail) {
    dt_ = dTnext;
    return fail;
  }

  /*
  // check solution integrity, fail if the solution has negative values
  for (int i = 0; i < soln_->SubVectors().size(); i++) {
    double min_var = 0;
    soln_->SubVectors()[i]->Data()->ViewComponent("cell")->MinValue(&min_var); 
    fail |=  min_var < -1e-12;
  }
  if (fail) {
    std::cout << "Step failed! Solution contains negative values. Restart! \n";
    *soln_->SubVector(0)->Data() = *S_->GetFieldData("phase1_pressure");
    *soln_->SubVector(1)->Data() = *S_->GetFieldData("phase2_saturation");
    return fail;
  }
  */

  nl_itrs_++;

  dt_ = dTnext;

  return fail;
}


int MPCoupled_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu)
{
  int ierr = 0;
  std::stringstream ss;
  ss << num_mat_;
  Amanzi::timer_manager.start("ApplyPreconditioner");
  if (cpr_enhanced_ || block_factorization_) {
    std::string file_name = "matrixA_" + ss.str() + ".txt";
    EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *comb_tree_op_->Op()->A());
    ierr = solver_comb_->ApplyInverse(*u, *pu);
    ln_itrs_ += solver_comb_->num_itrs();
  }
  else {
    std::string file_name = "matrixA_" + ss.str() + ".txt";
    EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());
    ierr = solver_tree_->ApplyInverse(*u, *pu);
    ln_itrs_ += solver_tree_->num_itrs();
  }
  Amanzi::timer_manager.stop("ApplyPreconditioner");
  std::string file_name = "rhs_p_" + ss.str() + ".txt";
  EpetraExt::MultiVectorToMatlabFile(file_name.c_str(), *u->SubVector(0)->Data()->ViewComponent("cell"));
  file_name = "rhs_s_" + ss.str() + ".txt";
  EpetraExt::MultiVectorToMatlabFile(file_name.c_str(), *u->SubVector(1)->Data()->ViewComponent("cell"));
  //std::cout << "rhs: \n";
  //u->Print(std::cout);
  //std::cout << "du: \n";
  //pu->Print(std::cout);
  //MPI_Comm mpi_comm(MPI_COMM_WORLD);
  //Amanzi::timer_manager.parSync(mpi_comm);
  //Amanzi::timer_manager.print();
  return ierr;
}

double MPCoupled_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) 
{
  double u1_inf, u2_inf;
  double du1_inf, du2_inf;
  const Epetra_MultiVector& du1_cell = *du->SubVector(0)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& du2_cell = *du->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& u1_cell = *u->SubVector(0)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& u2_cell = *u->SubVector(1)->Data()->ViewComponent("cell");

  du1_cell.NormInf(&du1_inf);
  du2_cell.NormInf(&du2_inf);
  u1_cell.NormInf(&u1_inf);
  u2_cell.NormInf(&u2_inf);
  double relative_du_inf = std::max(du1_inf/u1_inf, du2_inf/u2_inf);
  return relative_du_inf;

  //double l2_res;
  //du->Norm2(&l2_res);
  //return l2_res;

  /*
  double u_inf, du_inf;
  double du_l2;
  //u->NormInf(&u_inf);
  du->NormInf(&du_inf);
  du->Norm2(&du_l2);
  return du_l2;
  */
}


void MPCoupled_PK::CommitStep(double t_old, double t_new) {
  if (use_total_mobility_) {
    tot_phase_pk_->CommitStep(t_old, t_new);
  } else {
    phase1_pk_->CommitStep(t_old, t_new);
  }
  phase2_pk_->CommitStep(t_old, t_new);
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCoupled_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
           Teuchos::RCP<const TreeVector> u,
           Teuchos::RCP<TreeVector> du)
{
  Teuchos::RCP<CompositeVector> s_next = Teuchos::rcp(new CompositeVector(*u->SubVector(1)->Data()));
  s_next->Update(-1.0, *du->SubVector(1)->Data(), 1.0);
  ClipSaturation(s_next, 1e-8);
  du->SubVector(1)->Data()->Update(1.0, *u->SubVector(1)->Data(), -1.0, *s_next, 0.0);
}


void MPCoupled_PK::ClipSaturation(Teuchos::RCP<CompositeVector> s, double tol)
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

void MPCoupled_PK::ComputeAPinv()
{
  TreeVector unit_v(*soln_);
  TreeVector column_v(*soln_);
  int ncells = (*unit_v.SubVector(0)->Data()->ViewComponent("cell")).MyLength();
  //const char *file_name = "matrixA.txt";
  //EpetraExt::RowMatrixToMatlabFile(file_name, *tree_op_->A());
  std::stringstream ss;
  ss << num_mat_;
  std::string ln_str = ss.str();
  std::string file_out = "Pinv_" + ln_str;
  FILE* mat_out;
  mat_out = fopen(file_out.c_str(), "a");
  //mat_out.open("matrix_out.txt", ios::app);
  //mat_out << "Nonlinear iteration: " << nl_itrs_ << endl;
  //mat_out << scientific;
  //std::cout << std::scientific;
  //std::cout << "NL: " << nl_itrs_ << "; LN: " << ln_itrs_ << endl;
  //std::cout << *tree_op_->A() << "\n";
  for (int i = 0; i < 2*ncells; i++) {
    unit_v.PutScalar(0.0);
    column_v.PutScalar(0.0);
    if (i % 2 == 0) {
      Epetra_MultiVector& p_c = *unit_v.SubVector(0)->Data()->ViewComponent("cell");
      p_c[0][i/2] = 1.0;
      //std::cout << "modify pressure \n";
      //unit_v.Print(std::cout);
      if (cpr_enhanced_) {
        comb_tree_op_->ApplyInverse(unit_v, column_v);
      } else {
        tree_op_->ApplyInverse(unit_v, column_v);
      }
      //tree_op_->Apply(unit_v, column_v);
      Epetra_MultiVector& out_p_c = *column_v.SubVector(0)->Data()->ViewComponent("cell");
      Epetra_MultiVector& out_s_c = *column_v.SubVector(1)->Data()->ViewComponent("cell");
      //std::cout << "Printing column " << i << ": \n";
      for (int j = 0; j < ncells; j++){
        fprintf(mat_out, "%.8e %.8e ", out_p_c[0][j], out_s_c[0][j]);
        //std::cout << out_p_c[0][j] << " " << out_s_c[0][j] << " ";
        //mat_out << out_p_c[0][j] << " " << out_s_c[0][j] << " ";
      }
      fprintf(mat_out, "\n");
      //std::cout << endl;
      //mat_out << endl;
    }
    else {
      Epetra_MultiVector& s_c = *unit_v.SubVector(1)->Data()->ViewComponent("cell");
      s_c[0][(i-1)/2] = 1.0;
      //std::cout << "modify sat \n";
      //unit_v.Print(std::cout);
      tree_op_->ApplyInverse(unit_v, column_v);
      //tree_op_->Apply(unit_v, column_v);
      Epetra_MultiVector& out_p_c = *column_v.SubVector(0)->Data()->ViewComponent("cell");
      Epetra_MultiVector& out_s_c = *column_v.SubVector(1)->Data()->ViewComponent("cell");
      //std::cout << "Printing column " << i << ": \n";
      for (int j = 0; j < ncells; j++){
        fprintf(mat_out, "%.8e %.8e ", out_p_c[0][j], out_s_c[0][j]);
        //std::cout << out_p_c[0][j] << " " << out_s_c[0][j] << " ";
        //mat_out << out_p_c[0][j] << " " << out_s_c[0][j] << " ";
      }
      fprintf(mat_out, "\n");
      //std::cout << endl;
      //mat_out << endl;
    }
  }
  fflush(mat_out);
  fclose(mat_out);
  //mat_out.close();
}

}
}
