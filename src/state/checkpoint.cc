/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Checkpointing for state.

------------------------------------------------------------------------- */

#include "checkpoint.hh"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <iostream>
#include <iomanip>


namespace Amanzi {

Checkpoint::Checkpoint (Teuchos::ParameterList& plist, Epetra_MpiComm* comm) :
    IOEvent(plist, comm) {
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Checkpoint     ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_,this);

  // Set up the HDF5
  checkpoint_output_ = Teuchos::rcp(new HDF5_MPI(*comm));
  checkpoint_output_->setTrackXdmf(false);
}


// this constructor makes an object that will not create any output
Checkpoint::Checkpoint(): IOEvent() {}

// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Checkpoint::ReadParameters_() {
  filebasename_ = plist_.get<string>("file name base","checkpoint");
  filenamedigits_ = plist_.get<int>("file name digits", 5);
};

void Checkpoint::CreateFile(const int cycle) {
  // create the restart file
  std::stringstream oss;
  oss.flush();
  oss << filebasename_;
  oss.fill('0');
  oss.width(filenamedigits_);
  oss << std::right << cycle;
  checkpoint_output_->createDataFile(oss.str());
};


void Checkpoint::WriteVector(const Epetra_MultiVector& vec,
        const std::vector<std::string>& names ) const {
  if (names.size() < vec.NumVectors()) {
    Errors::Message m("Amanzi::Checkpoint::write_vector... not enough names were specified for the the components of the multi vector");
    Exceptions::amanzi_throw(m);
  }
  for (int i=0; i< vec.NumVectors(); i++) {
    checkpoint_output_->writeCellDataReal( *vec(i), names[i] );
  }
};

void Checkpoint::WriteAttributes(double time, double dt, int cycle) const {
  checkpoint_output_->writeAttrReal(time, "time");
  checkpoint_output_->writeAttrReal(dt, "dt");
  checkpoint_output_->writeAttrInt(cycle, "cycle");
  checkpoint_output_->writeAttrInt(comm_->NumProc(), "mpi_comm_world_rank");
};

} // namespace
