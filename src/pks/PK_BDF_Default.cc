/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon
*/

namespace Amanzi {

PK_BDF_Default::PK_BDF_Default(const Comm_ptr_type& comm,
        Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& glist,
        const Teuchos::RCP<State>& S)
  : PK_Default(comm, pk_tree, glist, S) {}


double
PK_BDF_Default::getDt()
{
  if (!strongly_coupled_)
    return S_->Get<double>("dt_internal", Tag(name_));
  else
    return -1.;
}


void
PK_BDF_Default::setDt(double dt)
{
  if (!strongly_coupled_) S_->Assign("dt_internal", Tag(name_), name_, dt);
}


void
PK_BDF_Default::parseParameterList()
{
  // preconditioner assembly
  assemble_preconditioner_ = plist_->get<bool>("assemble preconditioner", true);
  strongly_coupled_ = plist_->get<bool>("strongly coupled PK", false);
}


void
PK_BDF_Default::setup()
{
  if (!strongly_coupled_) {
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");

    // check if continuation method and require continuation parameter
    // -- ETC Note this needs fixed if more than one continuation method used
    if (bdf_plist.isSublist("continuation parameters")) {
      S_->Require<double>("continuation_parameter", Tag(name_), name_);
    }

    // require data for checkpointing timestep size
    S_->Require<double>("dt_internal", Tag(name_), name_);
  }
};


void
PK_BDF_Default::initialize()
{
  if (!strongly_coupled_) {
    // set up the timestepping algorithm
    // -- construct the time integrator
    //   Note, this is done here and not in setup because solution is not ready in setup
    auto bdf_plist = Teuchos::sublist(plist_, "time integrator");
    bdf_plist->sublist("verbose object")
      .setParametersNotAlreadySet(plist_->sublist("verbose object"));
    bdf_plist->sublist("verbose object").set("name", getName() + "_TI");

    solution_ = Teuchos::rcp(new TreeVector(getSolutionSpace()));
    State_to_Solution(tag_next_, *solution_);
    time_stepper_ =
      Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf_plist, solution_, S_));

    double dt_init = time_stepper_->initial_timestep();
    S_->Assign("dt_internal", Tag(name_), name_, dt_init);
    S_->GetRecordW("dt_internal", Tag(name_), name_).set_initialized();

    // -- initialize continuation parameter if needed.
    if (S_->HasRecord("continuation_parameter", Tag(name_))) {
      S_->Assign("continuation_parameter", Tag(name_), name_, (double)1.);
      S_->GetRecordW("continuation_parameter", Tag(name_), name_).set_initialized();
    }

    // -- initialize time derivative
    auto soln_dot = Teuchos::rcp(new TreeVector(*solution_));
    soln_dot->putScalar(0.0);

    // -- set initial state
    time_stepper_->SetInitialState(S_->get_time(), solution_, soln_dot);
  }
};


bool
PK_BDF_Default::advanceStep(double t_old, double t_new, bool reinit)
{
  double dt = t_new - t_old;
  Teuchos::OSTab out = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;
  State_to_Solution(Tags::NEXT, *solution_);

  // take a bdf timestep
  // Three dts:
  // --  dt is the requested timestep size.  It must be less than or equal to...
  // --  dt_internal is the max valid dt, and is set by physics/solvers
  // --  dt_solver is what the solver wants to do
  double dt_internal = S_->Get<double>("dt_internal", Tag(name_));

  // NOTE, still a bug in amanzi#685, despite fixes in amanzi#694, so this assertion still fails --ETC
  // AMANZI_ASSERT(dt <= dt_internal + 2.e-8); // roundoff

  double dt_solver = -1;
  bool fail = false;
  try {
    fail = time_stepper_->TimeStep(dt, dt_solver, solution_);

    if (!fail) {
      // check step validity
      bool valid = ValidStep();
      if (valid) {
        if (vo_->os_OK(Teuchos::VERB_LOW)) *vo_->os() << "successful advance" << std::endl;
        // update the timestep size
        if (dt_solver < dt_internal && dt_solver >= dt) {
          // We took a smaller step than we recommended, and it worked fine (not
          // suprisingly).  Likely this was due to constraints from other PKs or
          // vis.  Do not reduce our recommendation.
        } else {
          dt_internal = dt_solver;
        }
      } else {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "successful advance, but not valid" << std::endl;
        time_stepper_->CommitSolution(dt_internal, solution_, valid);
        dt_internal = 0.5 * dt_internal;
        // when including Valid here, make fail = true refs #110
      }
    } else {
      if (vo_->os_OK(Teuchos::VERB_LOW)) *vo_->os() << "unsuccessful advance" << std::endl;
      // take the decreased timestep size
      dt_internal = dt_solver;
    }

    S_->Assign("dt_internal", Tag(name_), name_, dt_internal);
  } catch (Errors::TimeStepCrash& e) {
    // inject more information into the crash message
    std::stringstream msg_str;
    msg_str << "TimeStepCrash in PK: \"" << getName() << "\"" << std::endl
            << "  at t = " << t_old << " with dt = " << dt << std::endl
            << "  error message: " << std::endl
            << std::endl
            << e.what() << std::endl
            << std::endl;
    Errors::TimeStepCrash msg(msg_str.str());
    Exceptions::amanzi_throw(msg);
  }
  return fail;
};


void
PK_BDF_Default::commitStep(double t_old, double t_new, const Tag& tag)
{
  if (tag == tag_next_) {
    double dt = t_new - t_old;
    if (time_stepper_ != Teuchos::null && dt > 0) {
      time_stepper_->CommitSolution(dt, solution_, true);
    }
  }
}


void
PK_BDF_Default::updateContinuationParameter(double lambda)
{
  S_->Assign("continuation_parameter", Tag(name_), name_, lambda);
  markChangedSolution();
}



} // namespace Amanzi
