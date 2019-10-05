#include "MPMC_PK.hh"
#include "EpetraExt_MultiVectorOut.h"

namespace Amanzi {
namespace Multiphase {

MPMC_PK::MPMC_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln):
                    soln_(soln), S_(S), pk_tree_(pk_tree), glist_(global_list),
                    passwd_("state")
{
  if (global_list->isSublist("Cycle Driver")) {
    ti_list_ = global_list->sublist("Cycle Driver").sublist("time integrator");   
  } else {
    ASSERT(0);
  }

  if (global_list->isSublist("MPMC Specs")) {
    mpmc_list_ = Teuchos::rcp(new Teuchos::ParameterList(global_list->sublist("MPMC Specs")));
  } else {
    ASSERT(0);
  }
  jacobian_type_ = mpmc_list_->get<std::string>("Jacobian type", "analytic");

  // We also need miscaleneous sublists
  if (global_list->isSublist("Preconditioners")) {
    preconditioner_list_ = global_list->sublist("Preconditioners");
  } else {
    Errors::Message msg("MPMC PK: input XML does not have <Preconditioners> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list->isSublist("Solvers")) {
    linear_operator_list_ = global_list->sublist("Solvers");
  } else {
    Errors::Message msg("MPMC PK: input XML does not have <Solvers> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  /*
  mesh_ = S->GetMesh();
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);

  p1_ = Teuchos::rcp(new CompositeVector(*cvs));
  s2_ = Teuchos::rcp(new CompositeVector(*cvs));
  p1_->PutScalar(1.0);
  s2_->PutScalar(0.8);
  */
  p1_tree_ = Teuchos::rcp(new TreeVector());
  s2_tree_ = Teuchos::rcp(new TreeVector());
  fuga1_tree_ = Teuchos::rcp(new TreeVector());
  fuga2_tree_ = Teuchos::rcp(new TreeVector());
  soln_->PushBack(p1_tree_);
  soln_->PushBack(s2_tree_);
  soln_->PushBack(fuga1_tree_);
  soln_->PushBack(fuga2_tree_);
  comp_w_pk_ = Teuchos::rcp(new Comp_PK(pk_tree.sublist("Component 1"), mpmc_list_, S, soln_));
  comp_h_pk_ = Teuchos::rcp(new Comp_PK(pk_tree.sublist("Component 2"), mpmc_list_, S, soln_));
  phase_constraint_pk_ = Teuchos::rcp(new Phase_Constraint_PK(*mpmc_list_, S));
  std::cout << "Done MPMC_PK Constructor\n";

  num_mat_ = 0;
}

MPMC_PK::~MPMC_PK() {
  // Do nothing for now
}

void MPMC_PK::Initialize() {
  rhs_ = Teuchos::rcp(new TreeVector());
  p1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("pressure_w",passwd_)));
  s2_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("saturation_n",passwd_)));
  fuga1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("fugacity1",passwd_)));
  fuga2_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("fugacity2",passwd_)));
  p1_tree_->SetData(p1_);
  s2_tree_->SetData(s2_);
  fuga1_tree_->SetData(fuga1_);
  fuga2_tree_->SetData(fuga2_);
  comp_w_pk_->Initialize();
  comp_h_pk_->Initialize();
  phase_constraint_pk_->Initialize();

  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  Teuchos::RCP<TreeVectorSpace> cvs_as_tv = Teuchos::rcp(new TreeVectorSpace(cvs));
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);
  tvs->PushBack(cvs_as_tv);

  tree_op_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  //tree_op_ = new Operators::TreeOperator(tvs);

  // Init time interval
  ProcessSublistTimeInterval(ti_list_, ti_specs_generic_);
 
  ti_specs_generic_.T0  = ti_list_.get<double>("start interval time", 0);
  ti_specs_generic_.dT0 = ti_list_.get<double>("initial time step", 1);

  double T0 = ti_specs_generic_.T0;
  double dT0 = ti_specs_generic_.dT0;

  dT = dT0;
  dTnext = dT0;

  std::cout<<"T0 "<<T0<<" dT0 "<<dT0<<"\n";

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
    //Teuchos::ParameterList bdf1_list = rp_list_.sublist(ti_method_name).sublist("BDF1");
    Teuchos::ParameterList bdf1_list = ti_specs_generic_.ti_list_ptr_->sublist("BDF1");
    //std::cout<<bdf1_list<<"\n";
    //if (! bdf1_list.isSublist("VerboseObject"))
    //    bdf1_list.sublist("VerboseObject") = mp_list_.sublist("VerboseObject");

    bdf1_dae = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  }
  std::cout << "MPMC_PK: done initialization\n";
}

bool MPMC_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  dT = t_new - t_old;
  dT_actual = dT;
  double time = S_->time();
  if (time >= 0.0) T_physics = time;
  time = T_physics;

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
      return fail;
    }
  }

  // tell the caller what time step we actually took
  dT_actual = dT;
  
  //dt_tuple times(time, dT);
  //ti_specs_->dT_history.push_back(times);

  ti_specs_->num_itrs++;

  dT = dTnext;

  return fail;
}

void MPMC_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  num_mat_++;
  std::cout << "Solution: \n";
  up->Print(std::cout);
  //std::cout << "MPMC_PK::UpdatePreconditioner\n";
  //comp_w_pk_->UpdatePreconditioner(t, up->SubVector(0), h);
  //comp_h_pk_->UpdatePreconditioner(t, up->SubVector(1), h);
  if (jacobian_type_ == "analytic")
  {
    comp_w_pk_->UpdatePreconditioner(t, up, h);
    comp_h_pk_->UpdatePreconditioner(t, up, h);

    //std::cout << *comp_w_pk_->op_pres_prec()->global_operator()->A() << "\n";
    //std::cout << *comp_w_pk_->op_sat_prec()->global_operator()->A() << "\n";
    //std::cout << *comp_h_pk_->op_pres_prec()->global_operator()->A() << "\n";
    //std::cout << *comp_h_pk_->op_sat_prec()->global_operator()->A() << "\n";
    tree_op_->SetOperatorBlock(0,0,comp_w_pk_->op_pres_prec()->global_operator());
    tree_op_->SetOperatorBlock(0,1,comp_w_pk_->op_sat_prec()->global_operator());
    tree_op_->SetOperatorBlock(0,2,comp_w_pk_->op_fuga_prec()->global_operator());
    tree_op_->SetOperatorBlock(1,0,comp_h_pk_->op_pres_prec()->global_operator());
    tree_op_->SetOperatorBlock(1,1,comp_h_pk_->op_sat_prec()->global_operator());
    tree_op_->SetOperatorBlock(1,3,comp_h_pk_->op_fuga_prec()->global_operator()); 
  } else if (jacobian_type_ == "numerical")
  {
    comp_w_pk_->NumericalJacobian(0.0, t, up, 1e-12);
    comp_h_pk_->NumericalJacobian(0.0, t, up, 1e-12);

    tree_op_->SetOperatorBlock(0,0,comp_w_pk_->Ops()[0]->global_operator());
    tree_op_->SetOperatorBlock(0,1,comp_w_pk_->Ops()[1]->global_operator());
    tree_op_->SetOperatorBlock(0,2,comp_w_pk_->Ops()[2]->global_operator());
    tree_op_->SetOperatorBlock(0,3,comp_w_pk_->Ops()[3]->global_operator());
    tree_op_->SetOperatorBlock(1,0,comp_h_pk_->Ops()[0]->global_operator());
    tree_op_->SetOperatorBlock(1,1,comp_h_pk_->Ops()[1]->global_operator());
    tree_op_->SetOperatorBlock(1,2,comp_h_pk_->Ops()[2]->global_operator());
    tree_op_->SetOperatorBlock(1,3,comp_h_pk_->Ops()[3]->global_operator());

  }

  phase_constraint_pk_->UpdatePreconditioner(t, up, h);
  tree_op_->SetOperatorBlock(2,0,phase_constraint_pk_->op_prec5()->global_operator());
  tree_op_->SetOperatorBlock(2,1,phase_constraint_pk_->op_p2_sat_prec()->global_operator());
  tree_op_->SetOperatorBlock(2,2,phase_constraint_pk_->op_prec3()->global_operator());
  tree_op_->SetOperatorBlock(2,3,phase_constraint_pk_->op_prec4()->global_operator()); 
  tree_op_->SetOperatorBlock(3,2,phase_constraint_pk_->op_prec1()->global_operator()); 
  tree_op_->SetOperatorBlock(3,3,phase_constraint_pk_->op_prec2()->global_operator()); 
  tree_op_->SetOperatorBlock(3,1,phase_constraint_pk_->op_p1_sat_prec()->global_operator()); 

  /*
  std::cout << *phase_constraint_pk_->op_prec5()->global_operator()->A() << "\n";
  std::cout << *phase_constraint_pk_->op_p2_sat_prec()->global_operator()->A() << "\n";
  std::cout << *phase_constraint_pk_->op_prec3()->global_operator()->A() << "\n";
  std::cout << *phase_constraint_pk_->op_prec4()->global_operator()->A() << "\n";
  std::cout << *phase_constraint_pk_->op_prec1()->global_operator()->A() << "\n";
  std::cout << *phase_constraint_pk_->op_prec2()->global_operator()->A() << "\n";
  std::cout << *phase_constraint_pk_->op_p1_sat_prec()->global_operator()->A() << "\n";
  */

  tree_op_->SymbolicAssembleMatrix();
  tree_op_->AssembleMatrix();
  
  std::stringstream ss;
  ss << num_mat_;
  std::string file_name = "matrixA_" + ss.str() + ".txt";
  EpetraExt::RowMatrixToMatlabFile(file_name.c_str(), *tree_op_->A());
  tree_op_->InitPreconditioner(ti_specs_->preconditioner_name, preconditioner_list_);
  //Teuchos::ParameterList tmp_list;
  //tree_op_->InitPreconditioner("identity", tmp_list);
  //std::cout << "after init prec: " << *tree_op_->A() << "\n";
}

void MPMC_PK::Functional(double t_old, double t_new, 
                         Teuchos::RCP<TreeVector> u_old,
                         Teuchos::RCP<TreeVector> u_new,
                         Teuchos::RCP<TreeVector> f) 
{
  comp_w_pk_->Functional(t_old, t_new, u_old, u_new, f->SubVector(0));
  comp_h_pk_->Functional(t_old, t_new, u_old, u_new, f->SubVector(1));
  phase_constraint_pk_->Functional(t_old, t_new, u_old, u_new, f);
  rhs_ = f;
  /*
  std::cout << "MPMC_PK::Functional\n";
  std::cout << *f->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  std::cout << *f->SubVector(1)->Data()->ViewComponent("cell") << "\n";
  std::cout << *f->SubVector(2)->Data()->ViewComponent("cell") << "\n";
  std::cout << *f->SubVector(3)->Data()->ViewComponent("cell") << "\n";
  */
}


int MPMC_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu)
{
  //std::cout << "MPMC_PK: ApplyPreconditioner\n";
  //std::cout << "RHS: \n";
  //std::cout << *u->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  //std::cout << *u->SubVector(1)->Data()->ViewComponent("cell") << "\n";  
  //std::cout << *u->SubVector(2)->Data()->ViewComponent("cell") << "\n";
  //std::cout << *u->SubVector(3)->Data()->ViewComponent("cell") << "\n";   

  AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::TreeOperator, TreeVector, TreeVectorSpace> >
      solver = factory.Create("GMRES", linear_operator_list_, tree_op_);
  solver->ApplyInverse(*u, *pu);

  //std::cout << "saturation correction: " << *pu->SubVector(1)->Data()->ViewComponent("cell") << "\n";

  //tree_op_->ApplyInverse(*u, *pu);
  //std::cout << "Solution: \n";
  //std::cout << *pu->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  //std::cout << *pu->SubVector(1)->Data()->ViewComponent("cell") << "\n";

  /*
  Teuchos::RCP<TreeVector> out = Teuchos::rcp(new TreeVector(*u));
  out->PutScalar(0.0);
  tree_op_->Apply(*pu, *out);
  std::cout << "apply to get functional back: \n";
  std::cout << *out->SubVector(0)->Data()->ViewComponent("cell") << "\n";
  std::cout << *out->SubVector(1)->Data()->ViewComponent("cell") << "\n";
  */
}

double MPMC_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) 
{
  //std::cout << "MPMC_PK: ErrorNorm\n";
  double u1_inf, u2_inf, u3_inf, u4_inf;
  double du1_inf, du2_inf, du3_inf, du4_inf;
  double du1_l2, du2_l2, du3_l2, du4_l2;
  const Epetra_MultiVector& du1_cell = *du->SubVector(0)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& du2_cell = *du->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& du3_cell = *du->SubVector(2)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& du4_cell = *du->SubVector(3)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& u1_cell = *u->SubVector(0)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& u2_cell = *u->SubVector(1)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& u3_cell = *u->SubVector(2)->Data()->ViewComponent("cell");
  const Epetra_MultiVector& u4_cell = *u->SubVector(3)->Data()->ViewComponent("cell");

  du1_cell.NormInf(&du1_inf);
  du2_cell.NormInf(&du2_inf);
  du3_cell.NormInf(&du3_inf);
  du4_cell.NormInf(&du4_inf);
  u1_cell.NormInf(&u1_inf);
  u2_cell.NormInf(&u2_inf);
  u3_cell.NormInf(&u3_inf);
  u4_cell.NormInf(&u4_inf);
  double relative_du_inf, du_inf, du_l2;
  relative_du_inf = std::max(du1_inf/u1_inf, du2_inf/u2_inf);
  relative_du_inf = std::max(relative_du_inf, du3_inf/u3_inf);
  relative_du_inf = std::max(relative_du_inf, du4_inf/u4_inf);
  du_inf = std::max(du1_inf, du2_inf);
  du_inf = std::max(du_inf, du3_inf);
  du_inf = std::max(du_inf, du4_inf);
  //std::cout << "Error l_inf: " << du_inf << "\n";

  du1_cell.Norm2(&du1_l2);
  du2_cell.Norm2(&du2_l2);
  du3_cell.Norm2(&du3_l2);
  du4_cell.Norm2(&du4_l2);

  du_l2 = std::max(du1_l2, du2_l2);
  du_l2 = std::max(du_l2, du3_l2);
  du_l2 = std::max(du_l2, du4_l2);

  //return du_inf;
  //return relative_du_inf;
  return du_l2;
}


void MPMC_PK::CommitStep(double t_old, double t_new) {
  comp_w_pk_->CommitStep(t_old, t_new);
  comp_h_pk_->CommitStep(t_old, t_new);
}


/* ******************************************************************
* Process generic time  interval sublist.
**************************************************************** */
void MPMC_PK::ProcessSublistTimeInterval(
  Teuchos::ParameterList& ti_list,  TI_Specs& ti_specs){

  //std::cout << "ti_list: \n" << ti_list << "\n";

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
void MPMC_PK::ProcessSublistTimeIntegration(
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
void MPMC_PK::ProcessStringTimeIntegration(const std::string name, int* method)
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
std::string MPMC_PK::FindStringPreconditioner(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("preconditioner")) {
    name = list.get<std::string>("preconditioner");
  } else {
    msg << "MPMC PK: parameter <preconditioner> is missing either in TI or LS list.";
    Exceptions::amanzi_throw(msg);
  }

  if (! preconditioner_list_.isSublist(name)) {
    msg << "MPMC PK: preconditioner \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Find string for the linear solver
**************************************************************** */
std::string MPMC_PK::FindStringLinearSolver(const Teuchos::ParameterList& list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("linear solver")) {
    name = list.get<std::string>("linear solver");
  } else {
    msg << "MPMC PK: time integrator does not define <linear solver>.";
    Exceptions::amanzi_throw(msg);
  }

  if (! linear_operator_list_.isSublist(name)) {
    msg << "MPMC PK: linear solver \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Process string for error control options
**************************************************************** */
void MPMC_PK::ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control)
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
