/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

IO event -- base class for reading or writing data.  Mostly just manages when
to do the I/O.

------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_IO_EVENT_HH_
#define AMANZI_STATE_IO_EVENT_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"

namespace Amanzi {

class TimeStepManager;

class IOEvent : public Teuchos::VerboseObject<IOEvent> {

 public:
  IOEvent(Teuchos::ParameterList& plist, Epetra_MpiComm *comm);
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

  // the Epetra communicator
  Epetra_MpiComm *comm_;

};

} // namespace

#endif
