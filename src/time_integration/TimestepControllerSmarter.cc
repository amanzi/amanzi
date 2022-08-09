/*
  Time Integration

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Slightly smarter timestep control based upon a history of previous timesteps.
*/

#include "errors.hh"
#include "dbc.hh"
#include "State.hh"
#include "TimestepControllerSmarter.hh"

namespace Amanzi {

TimestepControllerSmarter::TimestepControllerSmarter(const std::string& name,
        Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S) :
    plist_(plist),
    count_increased_before_increase_(0),
    successive_increases_(0),
    last_fail_(0),
    name_(name),
    S_(S)
{
  // allocate space for state -- done manually because Setup() has already been called
  if (S_ != Teuchos::null) {
    std::string varname = name_+"_increase_factor";
    S_->Require<double>(varname, Tags::DEFAULT, name_);
    S_->SetPtr<double>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new double(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    increase_factor_ = S_->GetPtrW<double>(varname, Tags::DEFAULT, name_);

    varname = name_+"_successive_increases";
    S_->Require<int>(varname, Tags::DEFAULT, name_);
    S_->SetPtr<int>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new int(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    successive_increases_ = S_->GetPtrW<int>(varname, Tags::DEFAULT, name_);

    varname = name_+"_last_fail";
    S_->Require<int>(varname, Tags::DEFAULT, name_);
    S_->SetPtr<int>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new int(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    last_fail_ = S_->GetPtrW<int>(varname, Tags::DEFAULT, name_);

    varname = name_+"_growth_wait_after_fail";
    S_->Require<int>(varname, Tags::DEFAULT, name_);
    S_->SetPtr<int>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new int(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    growth_wait_after_fail_ = S_->GetPtrW<int>(varname, Tags::DEFAULT, name_);
  } else {
    increase_factor_ = Teuchos::rcp(new double(0));
    successive_increases_ = Teuchos::rcp(new int(0));
    last_fail_ = Teuchos::rcp(new int(0));
    growth_wait_after_fail_ = Teuchos::rcp(new int(0));
  }

  max_its_ = plist_.get<int>("max iterations");
  min_its_ = plist_.get<int>("min iterations");
  AMANZI_ASSERT(max_its_ > min_its_);
  AMANZI_ASSERT(min_its_ >= 0);

  reduction_factor_ = plist_.get<double>("time step reduction factor");
  AMANZI_ASSERT(reduction_factor_ >= 0.0);
  AMANZI_ASSERT(reduction_factor_ <= 1.0);

  (*increase_factor_) = plist_.get<double>("time step increase factor");
  increase_factor0_ = (*increase_factor_);
  max_increase_factor_ = plist_.get<double>("max time step increase factor", 10.);
  AMANZI_ASSERT((*increase_factor_) >= 1.0);

  max_dt_ = plist_.get<double>("max time step");
  min_dt_ = plist_.get<double>("min time step");

  (*growth_wait_after_fail_) = plist_.get<int>("growth wait after fail");
  growth_wait_after_fail0_ = (*growth_wait_after_fail_);
  count_increased_before_increase_ = plist_.get<int>("count before increasing increase factor");

}


double
TimestepControllerSmarter::get_timestep(double dt, int iterations) {
  double dt_next(dt);

  // iterations < 0 implies failed timestep
  if (iterations < 0) {
    (*last_fail_) = 0;
    if ((*successive_increases_) > 0) {
      // Last time through we grew the timestep, and it failed.  Wait a bit
      // before growing the timestep.
      (*growth_wait_after_fail_)++;
    }

    (*increase_factor_) = std::min(increase_factor0_, max_increase_factor_);
    (*successive_increases_) = 0;

    dt_next = dt * reduction_factor_;

  } else {
    (*last_fail_)++;

    if ((*successive_increases_) > 0) {
      // Last time through we grew the timestep, and it was successful.  Reset
      // the wait counter.
      (*growth_wait_after_fail_) = growth_wait_after_fail0_;
    }

    if (iterations < min_its_) {
      if ((*last_fail_) > (*growth_wait_after_fail_)) { // grow the timestep
        (*successive_increases_)++;

        // increase faster than geometric if we are growing the timestep repeatedly
        if ((*successive_increases_) > count_increased_before_increase_) {
          (*increase_factor_) = std::min((*increase_factor_) * increase_factor0_, max_increase_factor_);
        }

        dt_next = dt * (*increase_factor_);
      }
    } else {
      // reset the counter
      (*successive_increases_) = 0;
      (*increase_factor_) = std::min(increase_factor0_, max_increase_factor_);

      if (iterations > max_its_) { // decrease the timestep
        dt_next = dt * reduction_factor_;
      }
    }
  }

  // check max step size
  if (dt_next > max_dt_) dt_next = max_dt_;

  // check min step size
  if (dt_next < min_dt_) {
    if (iterations < 0) {
      Errors::TimeStepCrash msg;
      msg << "Timestep failed: dT less than minimum (" << min_dt_ << ").";
      Exceptions::amanzi_throw(msg);
    } else {
      dt_next = min_dt_;
    }
  }

  // check that if we have failed, our step size has decreased.
  if (iterations < 0) {
    if (dt - dt_next < 1.e-10) {
      Errors::TimeStepCrash msg("Timestep failed: dT change is too small (check reduction_factor).");
      Exceptions::amanzi_throw(msg);
    }
  }

  return dt_next;
}

} // namespace Amanzi
