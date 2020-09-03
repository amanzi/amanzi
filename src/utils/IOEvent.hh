/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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

.. _io-event-spec:
.. admonition:: io-event-spec

   * `"cycles start period stop`" ``[Array(int)]`` **optional**

      The first entry is the start cycle, the second is the cycle
      period, and the third is the stop cycle or -1, in which case there
      is no stop cycle. A visualization dump is written at such
      cycles that satisfy cycle = start + n*period, for n=0,1,2,... and
      cycle < stop if stop != -1.0.

   * `"cycles start period stop 0`" ``[Array(int)]`` **optional** 

      If multiple cycles start period stop parameters are needed, then use these
      parameters.  If one with 0 is found, then one with 1 is looked for, etc,
      until the Nth one is not found.

   * `"cycles`" ``[Array(int)]``  **optional**
  
      An array of discrete cycles that at which a visualization dump is
      written.

   * `"times start period stop`" ``[Array(double)]`` **optional** 

      The first entry is the start time, the second is the time period,
      and the third is the stop time or -1, in which case there is no
      stop time. A visualization dump is written at such times that
      satisfy time = start + n*period, for n=0,1,2,... and time < stop
      if stop != -1.0.

   * `"times start period stop units`" ``[string]`` **s** 

      Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    
    * `"times start period stop 0`" ``[Array(double)]`` **optional**

      If multiple start period stop parameters are needed, then use this these
      parameters with N=0,1,2.  If one with 0 is found, then one with 1 is
      looked for, etc, until the Nth one is not found.

    * `"times start period stop 0 units`" ``[string]`` **s** 

      Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
      See above for continued integer listings.

    * `"times`" ``[Array(double)]`` **optional** 

      An array of discrete times that at which a visualization dump
      shall be written.

    * `"times units`" ``[string]`` **s** 

      Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    
 */

#ifndef AMANZI_STATE_IO_EVENT_HH_
#define AMANZI_STATE_IO_EVENT_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Units.hh"

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
  void ValidUnitOrThrow_(const std::string&);

  Teuchos::ParameterList plist_;

  Utils::Units units_;

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
