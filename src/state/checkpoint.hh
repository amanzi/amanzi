/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Checkpointing for state.

------------------------------------------------------------------------- */


#ifndef AMANZI_STATE_CHECKPOINT_HH_
#define AMANZI_STATE_CHECKPOINT_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_Comm.h"

#include "hdf5mpi_mesh.hh"
#include "io_event.hh"

namespace Amanzi {

class Checkpoint : public IOEvent {

 public:
  Checkpoint(Teuchos::ParameterList& plist, Epetra_MpiComm *comm);
  Checkpoint(); // this object will not create any output

  // public interface for coordinator clients
  void CreateFile(int cycle);
  void WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const;
  void WriteAttributes(double time, double dt, int cycle) const;

  void set_filebasename(std::string base) { filebasename_ = base; }

 protected:
  void ReadParameters_();

  std::string filebasename_;
  int filenamedigits_;
  int restart_cycle_;

  Teuchos::RCP<Amanzi::HDF5_MPI> checkpoint_output_;
};

}
#endif
