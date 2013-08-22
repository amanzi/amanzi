/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

IO event -- base class for reading or writing data.  Mostly just manages when
to do the I/O.

------------------------------------------------------------------------- */

#include "TimeStepManager.hh"

#include "io_event.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Standard constructor.
// -----------------------------------------------------------------------------
IOEvent::IOEvent(Teuchos::ParameterList& plist, Epetra_MpiComm *comm) :
    plist_(plist), disabled_(false), comm_(comm) {
  ReadParameters_();
};

// -----------------------------------------------------------------------------
// Constructor for a disabled Event.
// -----------------------------------------------------------------------------
IOEvent::IOEvent(): disabled_(true) {}


bool IOEvent::is_disabled() const { return disabled_; }
void IOEvent::disable(bool disabled) { disabled_ = disabled; }

// -----------------------------------------------------------------------------
// Place vis times in the manager.
// -----------------------------------------------------------------------------
void IOEvent::RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm) {
  if (times_.size() != 0) {
    tsm->RegisterTimeEvent(times_.toVector());
  }
  if (times_sps_.size() != 0) {
    for (Teuchos::Array<Teuchos::Array<double> >::const_iterator sps=times_sps_.begin();
         sps!=times_sps_.end(); ++sps) {
      tsm->RegisterTimeEvent((*sps)[0], (*sps)[1], (*sps)[2]);
    }
  }
}

// -----------------------------------------------------------------------------
// Does vis need to dump the state?
// -----------------------------------------------------------------------------
bool IOEvent::DumpRequested(int cycle, double time) const {
  return DumpRequested(cycle) || DumpRequested(time);
}

bool IOEvent::DumpRequested(int cycle) const {
  if (!is_disabled()) {
    if (cycles_.size() > 0) {
      for (Teuchos::Array<int>::const_iterator i=cycles_.begin(); i!=cycles_.end(); ++i) {
        if (cycle == *i) {
          return true;
        }
      }
    }
    if (cycles_sps_.size() != 0) {
      for (Teuchos::Array<Teuchos::Array<int> >::const_iterator sps=cycles_sps_.begin();
           sps!=cycles_sps_.end(); ++sps) {
        if ( (((*sps)[2]<0) || (cycle<=(*sps)[2])) && ((*sps)[1]>0) ) {
          if ((*sps)[0]<=cycle) {
            int cycle_loc = cycle - (*sps)[0];
            if (cycle_loc % (*sps)[1] == 0) return true;
          }
        }
      }
    }
  }
  // if none of the conditions apply we do not write a visualization dump
  return false;
}

bool IOEvent::DumpRequested(double time) const {
  if (!is_disabled()) {
    if (times_.size() > 0) {
      for (Teuchos::Array<double>::const_iterator i=times_.begin(); i!=times_.end(); ++i) {
        if (Amanzi::near_equal(*i,time)) return true;
      }
    }
    if (times_sps_.size() != 0) {
      for (Teuchos::Array<Teuchos::Array<double> >::const_iterator sps=times_sps_.begin();
           sps!=times_sps_.end(); ++sps) {
        if ( Amanzi::near_equal(time, (*sps)[0]) ) return true;
        if ((time > (*sps)[0] ) && ( ((*sps)[2] == -1.0) || ( time <= (*sps)[2] ))) {
          double n_periods = floor( (time - (*sps)[0])/(*sps)[1] );
          double next_time = (*sps)[0] + n_periods*(*sps)[1];
          if (Amanzi::near_equal(time, next_time)) return true;
        }
      }
    }
  }
  return false;
}


void IOEvent::ReadParameters_() {
  if (plist_.isParameter("cycles start period stop") ) {
    cycles_sps_.push_back(plist_.get<Teuchos::Array<int> >("cycles start period stop"));
  }

  bool done(false);
  int count(0);
  while (!done) {
    std::stringstream pname;
    pname << "cycles start period stop " << count;
    if (plist_.isParameter(pname.str())) {
      cycles_sps_.push_back(plist_.get<Teuchos::Array<int> >(pname.str()));
      count++;
    } else {
      done = true;
    }
  }

  if (plist_.isParameter("cycles")) {
    cycles_ = plist_.get<Teuchos::Array<int> >("cycles");
  }

  if (plist_.isParameter("times start period stop")) {
    times_sps_.push_back(plist_.get<Teuchos::Array<double> >("times start period stop"));
  }

  done = false;
  count = 0;
  while (!done) {
    std::stringstream pname;
    pname << "times start period stop " << count;
    if (plist_.isParameter(pname.str())) {
      times_sps_.push_back(plist_.get<Teuchos::Array<double> >(pname.str()));
      count++;
    } else {
      done = true;
    }
  }

  if (plist_.isParameter("times")) {
    times_ = plist_.get<Teuchos::Array<double> >("times");
  }
}


} // namespace
