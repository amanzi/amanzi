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
#include "IOEvent.hh"
#include "ObservationData.hh"

namespace Amanzi {

class Checkpoint : public IOEvent {

 public:
  Checkpoint(Teuchos::ParameterList& plist, Epetra_MpiComm *comm);
  Checkpoint(); // this object will not create any output

  // public interface for coordinator clients
  void CreateFile(int cycle);
  void WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const;
  void WriteAttributes(double time, double dt, int cycle, int pos) const;
  void WriteAttributes(double time, double dt, int cycle) const;
  void WriteAttributes(double time, int cycle) const;
  void WriteObservations(ObservationData* obs_data);
  void SetFinal(bool fnl) { final_ = fnl; }
  bool IsFinal() { return final_; }
  void Finalize();

  void set_filebasename(std::string base) { filebasename_ = base; }

 protected:
  void ReadParameters_();

  std::string filebasename_;
  int filenamedigits_;
  int restart_cycle_;
  bool final_;

  // the Epetra communicator
  const Epetra_MpiComm *comm_;
  
  Teuchos::RCP<Amanzi::HDF5_MPI> checkpoint_output_;
};

}
#endif
