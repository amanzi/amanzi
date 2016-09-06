/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Checkpointing for state.
------------------------------------------------------------------------- */
#include <iostream>
#include <iomanip>

#include "boost/filesystem.hpp"

#include "checkpoint.hh"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Amanzi {

Checkpoint::Checkpoint(Teuchos::ParameterList& plist, Epetra_MpiComm* comm) :
    IOEvent(plist),
    comm_(comm) {
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
  final_ = false;
}


// this constructor makes an object that will not create any output
Checkpoint::Checkpoint(): IOEvent() {}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Checkpoint::ReadParameters_() {
  filebasename_ = plist_.get<std::string>("file name base","checkpoint");
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
  checkpoint_output_->open_h5file();

  if (final_) { 
    std::string ch_file = oss.str() + ".h5";

    std::stringstream oss_final;
    oss_final << filebasename_<<"_final.h5";       
    std::string ch_final =  oss_final.str();

    if (boost::filesystem::is_regular_file(ch_final.data()))
      boost::filesystem::remove(ch_final.data());
    boost::filesystem::create_hard_link(ch_file.data(), ch_final.data());
  }
};


void Checkpoint::Finalize() {
  checkpoint_output_->close_h5file();
}


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


// -----------------------------------------------------------------------------
// Write simple attributes.
// -----------------------------------------------------------------------------
void Checkpoint::WriteAttributes(double time, double dt, int cycle, int position) const {
  checkpoint_output_->writeAttrReal(time, "time");
  checkpoint_output_->writeAttrReal(dt, "dt");
  checkpoint_output_->writeAttrInt(cycle, "cycle");
  checkpoint_output_->writeAttrInt(position, "position");
  checkpoint_output_->writeAttrInt(comm_->NumProc(), "mpi_comm_world_rank");
};


void Checkpoint::WriteAttributes(double time, double dt, int cycle) const {
  checkpoint_output_->writeAttrReal(time, "time");
  checkpoint_output_->writeAttrReal(dt, "dt");
  checkpoint_output_->writeAttrInt(cycle, "cycle");
  checkpoint_output_->writeAttrInt(comm_->NumProc(), "mpi_comm_world_rank");
};


void Checkpoint::WriteAttributes(double time, int cycle) const {
  checkpoint_output_->writeAttrReal(time, "time");
  checkpoint_output_->writeAttrInt(cycle, "cycle");
  checkpoint_output_->writeAttrInt(comm_->NumProc(), "mpi_comm_world_rank");
};


// -----------------------------------------------------------------------------
// Write observations
// -----------------------------------------------------------------------------
void Checkpoint::WriteObservations(ObservationData* obs_data)
{
  if (obs_data == NULL) return;

  std::vector<std::string> labels = obs_data->observationLabels();
  int nlabels = labels.size();

  if (nlabels > 0) {
    // save names of observations and their number
    int ndata(0);
    int *nobs; 
    char **tmp_labels;

    nobs = (int*)malloc(nlabels * sizeof(int));
    tmp_labels = (char**)malloc(nlabels * sizeof(char*));

    for (int i = 0; i < nlabels; i++) {
      tmp_labels[i] = (char*)malloc((labels[i].size() + 1) * sizeof(char));
      strcpy(tmp_labels[i], labels[i].c_str());
      nobs[i] = (*obs_data)[labels[i]].size();
      ndata += nobs[i];
    }

    checkpoint_output_->writeDataString(tmp_labels, nlabels, "obs_names");
    checkpoint_output_->writeAttrInt(nobs, nlabels, "obs_numbers");

    // save observation values
    double *tmp_data = (double*)malloc(2 * ndata * sizeof(double));
    int m(0);
    for (int i = 0; i < nlabels; ++i) {
      std::vector<ObservationData::DataQuadruple> tmp = (*obs_data)[labels[i]];
      for (int k = 0; k < tmp.size(); ++k) {
        tmp_data[m++] = tmp[k].time;
        tmp_data[m++] = tmp[k].value;
      }
    }
    checkpoint_output_->writeAttrReal(tmp_data, 2 * ndata, "obs_values");

    // release memory
    for (int i = 0; i < nlabels; i++) free(tmp_labels[i]);
    free(tmp_labels);
    free(nobs);
    free(tmp_data);
  }
}

}  // namespace Amanzi
