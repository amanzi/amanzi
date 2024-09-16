/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Slightly smarter timestep controller based upon a history of previous timesteps.
/*!

This is based on `Timestep Controller Standard`_, but also tries to be a bit
smarter to avoid repeated increase/decrease loops where the step size
decreases, converges in few iterations, increases, but then fails again.  It
also tries to grow the step geometrically to more quickly recover from tricky
nonlinearities.

.. _timestep-controller-smarter-spec:
.. admonition:: timestep-controller-smarter-spec

   * `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the
      previous step took more than this.
   * `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the
      previous step took less than this.
   * `"timestep reduction factor`" ``[double]`` :math:`f_{reduction}`, reduce
     the previous timestep by this multiple.
   * `"timestep increase factor`" ``[double]`` :math:`f_{increase}`, increase
     the previous timestep by this multiple.  Note that this can be modified
     geometrically in the case of repeated successful steps.
   * `"max timestep increase factor`" ``[double]`` **10.** The max
     :math:`f_{increase}` will ever get.
   * `"growth wait after fail`" ``[int]`` Wait at least this many timesteps
     before attempting to grow the timestep after a failed timestep.
   * `"count before increasing increase factor`" ``[int]`` Require this many
     successive increasions before multiplying :math:`f_{increase}` by itself.


*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "dbc.hh"
#include "State.hh"
#include "TimestepControllerRecoverable.hh"

namespace Amanzi {

template <class Vector>
class TimestepControllerSmarter : public TimestepControllerRecoverable<Vector> {
 public:
  TimestepControllerSmarter(const std::string& name,
                            Teuchos::ParameterList& plist,
                            const Teuchos::RCP<State>& S);

 protected:
  // single method for timestep control
  double getTimestep_(double dt, int iterations, bool valid) override;

 protected:
  int max_its_;
  int min_its_;

  double reduction_factor_;

  double increase_factor0_;
  double max_increase_factor_;
  int count_increased_before_increase_;
  int growth_wait_after_fail0_;

  // State variables stored in RCP to keep in state for checkpointing
  Teuchos::RCP<double> increase_factor_;
  Teuchos::RCP<int> successive_increases_;
  Teuchos::RCP<int> last_fail_;
  Teuchos::RCP<int> growth_wait_after_fail_;

  using TimestepControllerRecoverable<Vector>::S_;
  using TimestepControllerRecoverable<Vector>::name_;
};


template <class Vector>
TimestepControllerSmarter<Vector>::TimestepControllerSmarter(const std::string& name,
        Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S)
  : TimestepControllerRecoverable<Vector>(name, plist, S),
    count_increased_before_increase_(0),
    successive_increases_(0),
    last_fail_(0)
{
  // allocate space for state -- done manually because Setup() has already been called
  if (S_ != Teuchos::null) {
    std::string varname = name_ + "_increase_factor";
    S_->template Require<double>(varname, Tags::DEFAULT, name_);
    S_->template SetPtr<double>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new double(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_vis(false);
    increase_factor_ = S_->template GetPtrW<double>(varname, Tags::DEFAULT, name_);

    varname = name_ + "_successive_increases";
    S_->template Require<int>(varname, Tags::DEFAULT, name_);
    S_->template SetPtr<int>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new int(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_vis(false);
    successive_increases_ = S_->template GetPtrW<int>(varname, Tags::DEFAULT, name_);

    varname = name_ + "_last_fail";
    S_->template Require<int>(varname, Tags::DEFAULT, name_);
    S_->template SetPtr<int>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new int(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_vis(false);
    last_fail_ = S_->template GetPtrW<int>(varname, Tags::DEFAULT, name_);

    varname = name_ + "_growth_wait_after_fail";
    S_->template Require<int>(varname, Tags::DEFAULT, name_);
    S_->template SetPtr<int>(varname, Tags::DEFAULT, name_, Teuchos::rcp(new int(0)));
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_checkpoint();
    S_->GetRecordW(varname, Tags::DEFAULT, name_).set_io_vis(false);
    growth_wait_after_fail_ = S_->template GetPtrW<int>(varname, Tags::DEFAULT, name_);
  } else {
    increase_factor_ = Teuchos::rcp(new double(0));
    successive_increases_ = Teuchos::rcp(new int(0));
    last_fail_ = Teuchos::rcp(new int(0));
    growth_wait_after_fail_ = Teuchos::rcp(new int(0));
  }

  max_its_ = plist.get<int>("max iterations");
  min_its_ = plist.get<int>("min iterations");
  AMANZI_ASSERT(max_its_ > min_its_);
  AMANZI_ASSERT(min_its_ >= 0);

  reduction_factor_ = plist.get<double>("timestep reduction factor");
  AMANZI_ASSERT(reduction_factor_ >= 0.0);
  AMANZI_ASSERT(reduction_factor_ <= 1.0);

  (*increase_factor_) = plist.get<double>("timestep increase factor");
  increase_factor0_ = (*increase_factor_);
  max_increase_factor_ = plist.get<double>("max timestep increase factor", 10.);
  AMANZI_ASSERT((*increase_factor_) >= 1.0);

  (*growth_wait_after_fail_) = plist.get<int>("growth wait after fail");
  growth_wait_after_fail0_ = (*growth_wait_after_fail_);
  count_increased_before_increase_ = plist.get<int>("count before increasing increase factor");
}


template <class Vector>
double
TimestepControllerSmarter<Vector>::getTimestep_(double dt, int iterations, bool valid)
{
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
          (*increase_factor_) =
            std::min((*increase_factor_) * increase_factor0_, max_increase_factor_);
        }

        dt_next = dt * (*increase_factor_);
      }
    } else {
      // reset the counter
      (*successive_increases_) = 0;
      (*increase_factor_) = std::min(increase_factor0_, max_increase_factor_);

      if (iterations > max_its_ || !valid) { // decrease the timestep
        dt_next = dt * reduction_factor_;
      }
    }
  }
  return dt_next;
}

} // namespace Amanzi


