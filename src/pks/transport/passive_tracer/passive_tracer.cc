/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   ------------------------------------------------------------------------- */
#include "advection_factory.hh"
#include "bdf2_time_integrator.hh"
#include "passive_tracer.hh"
#include "tabular-function.hh"

namespace Amanzi {
namespace Transport {

PassiveTracer::PassiveTracer(Teuchos::ParameterList& transport_plist,
                             Teuchos::ParameterList& FElist,
                             const Teuchos::RCP<TreeVector>& solution) :
  transport_plist_(transport_plist) {

  // require fields for the state and solution
  S->RequireField("concentration", "transport", AmanziMesh::CELL);
  S->GetField("concentration","transport")->set_io_vis(true);
  Teuchos::RCP<CompositeVector> conc = S->GetFieldData("concentration", "transport");
  solution->SetData(conc);
  solution_ = solution;

  // independent variables (not owned by this pk)
  S->RequireField("porosity", AmanziMesh::CELL, 1, true);
  S->RequireField("saturation_liquid", AmanziMesh::CELL, 1, true);
  S->RequireField("mass_flux", AmanziMesh::FACE, 1, true);
  S->RequireField("cell_volume", AmanziMesh::CELL, 1, true);

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = transport_plist_.sublist("advection");
  advection_ = advection_factory.create(advect_plist, S->Mesh());

  // bcs
  bcs_ = Teuchos::rcp(new std::vector< Teuchos::RCP<BoundaryFunction> >);
};

// initialize ICs
void PassiveTracer::Initialize(const Teuchos::RCP<State>& S) {
  process_parameter_list(S);

  // constant initial concentration
  C_ = transport_plist_.get<double>("Constant concentration", 1.0);
  S->GetFieldData("concentration", "transport")->PutScalar(C_);
  S->GetField("concentration", "transport")->set_initialized();

  // initialize the advection method
  advection_->set_num_dofs(1);

  state_to_solution(S, solution_);
  atol_ = transport_plist_.get<double>("Absolute error tolerance",1.0);
  rtol_ = transport_plist_.get<double>("Relative error tolerance",1.0);

  if (!transport_plist_.get<bool>("Strongly Coupled PK", false)) {
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(transport_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};

void PassiveTracer::process_parameter_list(const Teuchos::RCP<State>& S) {
  // global transport parameters
  cfl_ = transport_plist_.get<double>("CFL", 1.);

  // stupid hack, assumes max(v) ~ 1.
  dt_ = .1*cfl_;

  // extract list of lists of boundary conditions
  Teuchos::ParameterList BCs_list;
  BCs_list = transport_plist_.get<Teuchos::ParameterList>("Transport BCs");

  // populate the list of boundary influx functions
  int nBCs = BCs_list.get<int>("number of BCs");
  bcs_->clear();
  bcs_->reserve(nBCs);

  for (int n=0; n != nBCs; ++n) {
    char bc_char_name[10];
    sprintf(bc_char_name, "BC %d", n);
    string bc_name(bc_char_name);

    if (!BCs_list.isSublist(bc_name)) {
      Errors::Message msg;
      msg << "Boundary condition with name " << bc_char_name << " does not exist" << "\n";
      Exceptions::amanzi_throw(msg);
    }
    Teuchos::ParameterList BC_list = BCs_list.sublist(bc_name);  // A single sublist.

    int i = 0;
    char tcc_char_name[20];
    sprintf(tcc_char_name, "Component %d", i);
    string tcc_name(tcc_char_name);

    bool flag_BCX = false;
    if (BC_list.isParameter(tcc_name)) {
      flag_BCX = true;
      std::vector<std::string> regions, functions;
      std::vector<double> times, values;

      regions = BC_list.get<Teuchos::Array<std::string> >("Regions").toVector();
      times = BC_list.get<Teuchos::Array<double> >("Times").toVector();
      values = BC_list.get<Teuchos::Array<double> >(tcc_name).toVector();
      functions = BC_list.get<Teuchos::Array<std::string> >("Time Functions").toVector();

      int nfunctions = functions.size();  // convert strings to forms
      std::vector<TabularFunction::Form> forms(functions.size());
      for (int k=0; k != nfunctions; ++k) {
        forms[k] = (functions[k] == "Constant") ? TabularFunction::CONSTANT : TabularFunction::LINEAR;
      }

      Teuchos::RCP<Function> func;
      func = Teuchos::rcp(new TabularFunction(times, values, forms));
      Teuchos::RCP<BoundaryFunction> bnd_fun = Teuchos::rcp(new BoundaryFunction(S->Mesh()));
      bnd_fun->Define(regions, func);
      bcs_->push_back(bnd_fun);
    }
    if (!flag_BCX) {
      Errors::Message msg;
      msg << "Sublist BC X was not found.\n";
      Exceptions::amanzi_throw(msg);
    }
  }
  //print_statistics();
}



  
// Pointer copy of state to solution
void PassiveTracer::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) {
  soln->SetData(S->GetFieldData("concentration", "transport"));
};

// Pointer copy concentration fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void PassiveTracer::solution_to_state(const Teuchos::RCP<TreeVector>& soln,
        const Teuchos::RCP<State>& S) {
  S->SetData("concentration", "transport", soln->Data());
};

// -- advance using the BDF integrator
bool PassiveTracer::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);
  time_stepper_->commit_solution(h, solution_);
  return false;
};


} // namespace
} // namespace
