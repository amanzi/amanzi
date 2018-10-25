/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt

Checkpointing for state.

NOTE: Should make this class RAII
------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_CHECKPOINT_HH_
#define AMANZI_STATE_CHECKPOINT_HH_

#include "AmanziTypes.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "HDF5_MPI.hh"
#include "IOEvent.hh"

namespace Amanzi {

class Checkpoint : public IOEvent {

public:
  // standard output constructor
  Checkpoint(Teuchos::ParameterList &plist, Comm_ptr_type comm);

  // standard input constructor
  Checkpoint(const std::string &filename, Comm_ptr_type comm);

  // this object will not create any output
  Checkpoint();

  // public interface for coordinator clients
  void CreateFile(int cycle);

  template <typename T> void Write(const std::string &name, const T &t) const;

  template <typename T> void Read(const std::string &name, T &t) const {}

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
  Comm_ptr_type comm_;
  Teuchos::RCP<Amanzi::HDF5_MPI> checkpoint_output_;
};

template <>
inline void Checkpoint::Write<Epetra_Vector>(const std::string &name,
                                             const Epetra_Vector &t) const {
  checkpoint_output_->writeCellDataReal(t, name);
}

template <>
inline void Checkpoint::Write<double>(const std::string &name,
                                      const double &t) const {
  checkpoint_output_->writeAttrReal(t, name);
}

template <>
inline void Checkpoint::Write<int>(const std::string &name,
                                   const int &t) const {
  checkpoint_output_->writeAttrInt(t, name);
}

template <>
inline void Checkpoint::Read<Epetra_Vector>(const std::string &name,
                                            Epetra_Vector &t) const {
  checkpoint_output_->readData(t, name);
}

template <>
inline void Checkpoint::Read<double>(const std::string &name, double &t) const {
  checkpoint_output_->readAttrReal(t, name);
}

template <>
inline void Checkpoint::Read<int>(const std::string &name, int &t) const {
  checkpoint_output_->readAttrInt(t, name);
}

} // namespace Amanzi

#endif
