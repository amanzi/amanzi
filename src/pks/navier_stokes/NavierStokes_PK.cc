/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "PK_DomainFunctionFactory.hh"

// Amanzi::NavierStokes
#include "NavierStokes_PK.hh"

namespace Amanzi {
namespace NavierStokes {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
NavierStokes_PK::NavierStokes_PK(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln),
    passwd_("navier stokes")
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ns_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ns_list_, "time integrator", true);
   
  // domain name
  domain_ = ns_list_->get<std::string>("domain name", "domain");
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
NavierStokes_PK::NavierStokes_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const std::string& pk_list_name,
                                 Teuchos::RCP<State> S,
                                 const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    soln_(soln),
    passwd_("navier stokes")
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ns_list_ = Teuchos::sublist(pk_list, "Navier Stokes", true);
 
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ns_list_, "time integrator");

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void NavierStokes_PK::Setup(const Teuchos::Ptr<State>& S)
{
  dt_ = 0.0;
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  // primary fields
  // -- pressure
  if (!S->HasData("pressure")) {
    S->Require<CV_t, CVS_t>("pressure", Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist("pressure");
    elist.set<std::string>("evaluator name", "pressure");
    pressure_eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S->SetEvaluator("pressure", pressure_eval_);
  }

  // -- velocity
  std::vector<std::string> names = {"node", "face"};
  std::vector<AmanziMesh::Entity_kind> locations = {AmanziMesh::NODE, AmanziMesh::FACE};
  std::vector<int> ndofs = {dim, 1};

  if (!S->HasData("fluid_velocity")) {
    S->Require<CV_t, CVS_t>("fluid_velocity", Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist("fluid_velocity");
    elist.set<std::string>("evaluator name", "fluid_velocity");
    fluid_velocity_eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S->SetEvaluator("fluid_velocity", fluid_velocity_eval_);
  }

  // -- viscosity: if not requested by any PK, we request its constant value.
  if (!S->HasData("const_fluid_viscosity")) {
    S->Require<double>("const_fluid_viscosity", Tags::DEFAULT, "state");
  }
}


/* ******************************************************************
* This is a long but simple routine. It goes through flow parameter
* list and initializes various objects including those created during 
* the setup step.
****************************************************************** */
void NavierStokes_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Initialize miscalleneous defaults.
  // -- times
  double t_ini = S->get_time(); 
  dt_desirable_ = dt_;
  dt_next_ = dt_;

  // -- mesh dimensions
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // Create verbosity object to print out initialiation statistics.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ns_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("NavierStokes", vlist)); 

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os()<< "\nPK initialization started...\n";
  }

  // Create pointers to the primary flow field pressure.
  Teuchos::RCP<TreeVector> tmp_u = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> tmp_p = Teuchos::rcp(new TreeVector());
  soln_->PushBack(tmp_u);
  soln_->PushBack(tmp_p);
 
  soln_p_ = S->GetPtrW<CV_t>("pressure", Tags::DEFAULT, passwd_);
  soln_u_ = S->GetPtrW<CV_t>("fluid_velocity", Tags::DEFAULT, passwd_);
  soln_->SubVector(0)->SetData(soln_u_); 
  soln_->SubVector(1)->SetData(soln_p_); 

  // Initialize time integrator.
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  AMANZI_ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (! bdf1_list.isSublist("verbose object"))
      bdf1_list.sublist("verbose object") = ns_list_->sublist("verbose object");

  bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));

  // Initialize matrix and preconditioner
  // -- create elastic block
  Teuchos::ParameterList& tmp1 = ns_list_->sublist("operators")
                                          .sublist("elasticity operator");
  op_matrix_elas_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));
  op_preconditioner_elas_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));

  // -- create divergence block
  Teuchos::ParameterList& tmp2 = ns_list_->sublist("operators")
                                          .sublist("divergence operator");
  op_matrix_div_ = Teuchos::rcp(new Operators::PDE_Abstract(tmp2, mesh_));

  // -- create gradient block (transpose of divergence block)
  Teuchos::ParameterList& tmp3 = ns_list_->sublist("operators")
                                          .sublist("gradient operator");
  op_matrix_grad_ = Teuchos::rcp(new Operators::PDE_Abstract(tmp3, mesh_));

  // -- create accumulation term (velocity block, only nodal unknowns)
  Operators::Schema schema(AmanziMesh::NODE, 2);  // FIXME
  op_matrix_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(schema, op_matrix_elas_->global_operator()));
  op_preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(schema, op_preconditioner_elas_->global_operator()));

  // -- create convection term
  Teuchos::ParameterList& tmp4 = ns_list_->sublist("operators")
                                          .sublist("advection operator");
  op_matrix_conv_ = Teuchos::rcp(new Operators::PDE_Abstract(tmp4, op_matrix_elas_->global_operator()));
  op_preconditioner_conv_ = Teuchos::rcp(new Operators::PDE_Abstract(tmp4, op_preconditioner_elas_->global_operator()));

  // -- create pressure block (for preconditioner)
  op_mass_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));

  // -- matrix and preconditioner
  op_matrix_ = Teuchos::rcp(new Operators::TreeOperator(Teuchos::rcpFromRef(soln_->Map())));
  op_matrix_->set_operator_block(0, 0, op_matrix_elas_->global_operator());
  op_matrix_->set_operator_block(1, 0, op_matrix_div_->global_operator());
  op_matrix_->set_operator_block(0, 1, op_matrix_grad_->global_operator());

  op_preconditioner_ = Teuchos::rcp(new Operators::TreeOperator(Teuchos::rcpFromRef(soln_->Map())));
  op_preconditioner_->set_operator_block(0, 0, op_preconditioner_elas_->global_operator());
  op_preconditioner_->set_operator_block(1, 1, op_mass_->global_operator());

  // Create BC objects
  Teuchos::RCP<NavierStokesBoundaryFunction> bc;
  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(ns_list_->sublist("boundary conditions", true)));

  bcs_.clear();

  // -- velocity
  if (bc_list->isSublist("velocity")) {
    PK_DomainFunctionFactory<NavierStokesBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("velocity");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        auto schema_pc = op_preconditioner_elas_->global_schema_col();
        for (auto jt = schema_pc.begin(); jt != schema_pc.end(); ++jt) {
          AmanziMesh::Entity_kind kind;
          WhetStone::DOF_Type type;
          std::tie(kind, type, std::ignore) = *jt;

          bc = bc_factory.Create(spec, "no slip", kind, Teuchos::null);
          bc->set_bc_name("no slip");
          bc->set_type(type);
          bcs_.push_back(bc);
        }
      }
    }
  }

  // Populate matrix and preconditioner
  // -- setup phase
  double mu = S_->Get<double>("const_fluid_viscosity");
  op_matrix_elas_->global_operator()->Init();
  op_matrix_elas_->SetTensorCoefficient(mu);

  op_preconditioner_elas_->global_operator()->Init();
  op_preconditioner_elas_->SetTensorCoefficient(mu);

  // -- add boundary conditions for velocity
  op_bcs_.clear();
  schema = op_matrix_elas_->schema_col();
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    AmanziMesh::Entity_kind kind;
    WhetStone::DOF_Type type;
    std::tie(kind, type, std::ignore) = *it;

    auto bcx = Teuchos::rcp(new Operators::BCs(mesh_, kind, type));
    op_matrix_elas_->AddBCs(bcx, bcx);
    op_matrix_conv_->AddBCs(bcx, bcx);
    op_matrix_acc_->AddBCs(bcx, bcx);
    op_matrix_div_->AddBCs(bcx, bcx);
    op_matrix_grad_->AddBCs(bcx, bcx);

    op_preconditioner_elas_->AddBCs(bcx, bcx);
    op_preconditioner_conv_->AddBCs(bcx, bcx);
    op_preconditioner_acc_->AddBCs(bcx, bcx);

    op_bcs_.push_back(bcx);
  }

  // -- iniailize boundary conditions (memory allocation)
  UpdateSourceBoundaryData_(t_ini, t_ini);
  op_matrix_elas_->EnforceBCs(*soln_u_);

  // -- assemble phase
  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  op_preconditioner_elas_->UpdateMatrices();
  op_preconditioner_elas_->ApplyBCs(true, true, true);

  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  op_preconditioner_elas_->global_operator()->set_inverse_parameters(pc_name, *preconditioner_list_);

  CompositeVector vol(op_mass_->global_operator()->DomainMap());
  vol.PutScalar(1.0 / mu);
  op_mass_->AddAccumulationTerm(vol, 1.0, "cell");
  op_mass_->global_operator()->set_inverse_parameters("Diagonal", *preconditioner_list_);

  // -- generic linear solver for most cases
  solver_name_ = ti_list_->get<std::string>("linear solver");

  Teuchos::ParameterList inv_list;
  inv_list.set("preconditioning method", "block diagonal");
  if (ti_list_->isParameter("preconditioner enhancement")) {
    std::string tmp_solver = ti_list_->get<std::string>("preconditioner enhancement");
    inv_list.setParameters(linear_solver_list_->sublist(tmp_solver));
  }
  op_preconditioner_->set_inverse_parameters(inv_list);

  // summary of initialization
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << " TI:\"" << ti_method_name.c_str() << "\"" << std::endl
               << "matrix:\n" << op_matrix_->PrintDiagnostics() << std::endl
               << "precon:\n" << op_preconditioner_->PrintDiagnostics() << std::endl;

    // *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    // *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl << std::endl;

    *vo_->os() << vo_->color("green") << "Initialization of PK is complete, T=" 
               << units_.OutputTime(S_->get_time()) << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************* 
* Performs one time step from time t_old to time t_new either for
* steady-state or transient simulation.
******************************************************************* */
bool NavierStokes_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // save a copy of primary and conservative fields
  CompositeVector pressure_copy(S_->Get<CV_t>("pressure"));
  CompositeVector fluid_velocity_copy(S_->Get<CV_t>("fluid_velocity"));

  // initialization
  if (num_itrs_ == 0) {
    Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae_->SetInitialState(t_old, soln_, udot);

    UpdatePreconditioner(t_old, soln_, dt_);
    num_itrs_++;
  }

  // trying to make a step
  bool failed(false);
  failed = bdf1_dae_->TimeStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;

    // recover the original primary solution
    S_->GetW<CV_t>("pressure", Tags::DEFAULT, passwd_) = pressure_copy;
    pressure_eval_->SetChanged();

    S_->GetW<CV_t>("fluid_velocity", Tags::DEFAULT, passwd_) = fluid_velocity_copy;
    fluid_velocity_eval_->SetChanged();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Reverted pressure, fluid_velocity" << std::endl;

    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae_->CommitSolution(dt_, soln_);
  pressure_eval_->SetChanged();
  fluid_velocity_eval_->SetChanged();

  num_itrs_++;
  dt_ = dt_next_;
  
  return failed;
}
 

/* ******************************************************************* 
* Performs one time step from time t_old to time t_new either for
* steady-state or transient simulation.
******************************************************************* */
void NavierStokes_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  double tmp1, tmp2;
  soln_u_->Norm2(&tmp1);
  soln_p_->Norm2(&tmp2);
  *vo_->os() << "solution norms=" << tmp1 << " " << tmp2 << std::endl;
}

}  // namespace NavierStokes
}  // namespace Amanzi
