/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! IOEvent: base time/timestep control determing when in time to do something.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)

*/

/*!

The IOEvent is used for multiple objects that need to indicate simulation times or cycles on which to do something.

* `"cycles start period stop`" ``[Array(int)]`` 

    The first entry is the start cycle, the second is the cycle
    period, and the third is the stop cycle or -1, in which case there
    is no stop cycle. A visualization dump is written at such
    cycles that satisfy cycle = start + n*period, for n=0,1,2,... and
    cycle < stop if stop != -1.0.

* `"cycles start period stop N`" ``[Array(int)]`` 

    If multiple cycles start period stop parameters are needed, then
    use these parameters with N=0,1,2,...

* `"cycles`" ``[Array(int)]`` 
  
    An array of discrete cycles that at which a visualization dump is
    written.

* `"times start period stop`" ``[Array(double)]`` 

    The first entry is the start time, the second is the time period,
    and the third is the stop time or -1, in which case there is no
    stop time. A visualization dump is written at such times that
    satisfy time = start + n*period, for n=0,1,2,... and time < stop
    if stop != -1.0.  Note that all times units are in seconds.

* `"times start period stop n`" ``[Array(double)]``

    If multiple start period stop parameters are needed, then use this
    these parameters with n=0,1,2,..., and not the single `"times
    start period stop`" parameter.  Note that all times units are in
    seconds.

* `"times`" ``[Array(double)]`` 

    An array of discrete times that at which a visualization dump
    shall be written.  Note that all times units are in seconds.
 */

#ifndef AMANZI_STATE_IO_EVENT_HH_
#define AMANZI_STATE_IO_EVENT_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Amanzi {

class TimeStepManager;

class IOEvent : public Teuchos::VerboseObject<IOEvent> {
 public:
  IOEvent(Teuchos::ParameterList& plist);
  IOEvent(); // created with this constructor this object will not create any output

  void disable(bool disabled=true);
  bool is_disabled() const;

  // public interface for coordinator clients
  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);
  bool DumpRequested(int cycle, double time) const;
  bool DumpRequested(int cycle) const;
  bool DumpRequested(double time) const;

 protected:
  void ReadParameters_();

  Teuchos::ParameterList plist_;

  // Time step control -- when to do this i/o?
  Teuchos::Array<int> cycles_;
  Teuchos::Array<Teuchos::Array<int> > cycles_sps_;
  Teuchos::Array<double> times_;
  Teuchos::Array<Teuchos::Array<double> > times_sps_;

  // disable visualization dumps alltogether
  bool disabled_;
};

} // namespace Amanzi

#endif
