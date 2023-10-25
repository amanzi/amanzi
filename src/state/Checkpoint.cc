/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/* -------------------------------------------------------------------------
Amanzi

License:

Checkpointing for state.
------------------------------------------------------------------------- */
#include <filesystem>
#include <iostream>
#include <iomanip>

#include "Epetra_MpiComm.h"
#include "mpi.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Checkpoint.hh"
#include "State.hh"

namespace Amanzi {

Checkpoint::Checkpoint(bool single_file) : single_file_(single_file)
{
  output_["domain"] = Teuchos::rcp(new HDF5_MPI(Amanzi::getDefaultComm()));
  output_["domain"]->setTrackXdmf(false);
}

Checkpoint::Checkpoint(Teuchos::ParameterList& plist, const State& S)
  : IOEvent(plist), single_file_(true)
{
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Checkpoint     ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_, this);

  // Set up the HDF5 objects
  if (single_file_) {
    // Single_File style is one checkpoint file on "domain"'s comm
    Comm_ptr_type comm;
    if (S.HasMesh("domain")) {
      comm = S.GetMesh()->getComm();
    } else {
      comm = Amanzi::getDefaultComm();
    }
    output_["domain"] = Teuchos::rcp(new HDF5_MPI(comm));
    output_["domain"]->setTrackXdmf(false);

    // This requires that all meshes and data are defined on MPI_COMM_WORLD
    // Check to confirm!
    for (auto entry = S.mesh_begin(); entry != S.mesh_end(); ++entry) {
      const std::string& domain = entry->first;
      const auto& mesh = S.GetMesh(domain);
      if (!sameComm(*comm, *mesh->getComm())) {
        std::stringstream msg;
        msg << "Checkpointing: cannot use single file checkpointing when not all meshes are on "
               "MPI_COMM_WORLD (mesh \""
            << domain
            << "\" not on MPI_COMM_WORLD).  Using multi-file checkpointing.  (Hide this warning by "
               "setting \"single file checkpoint\" to \"false\" in the \"checkpoint\" list.";
        std::cerr << "WARNING: " << msg.str() << std::endl;
        single_file_ = false;
        break;
      }
    }
  }

  // NOTE: do not make this an 'else' clause!
  if (!single_file_) {
    for (auto domain = S.mesh_begin(); domain != S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);
      output_[domain->first] = Teuchos::rcp(new HDF5_MPI(mesh->getComm()));
      output_[domain->first]->setTrackXdmf(false);
    }
  }
}


// Constructor for reading
Checkpoint::Checkpoint(const std::string& file_or_dirname, const State& S) : IOEvent()
{
  // if provided a directory, use new style
  if (std::filesystem::is_directory(file_or_dirname)) {
    single_file_ = false;
  } else if (std::filesystem::is_regular_file(file_or_dirname)) {
    single_file_ = true;
  } else {
    Errors::Message message;
    message << "Checkpoint::Read: location \"" << file_or_dirname << "\" does not exist.";
    Exceptions::amanzi_throw(message);
  }

  // create the readers
  if (single_file_) {
    auto comm = S.GetMesh()->getComm();
    output_["domain"] = Teuchos::rcp(new HDF5_MPI(comm, file_or_dirname));
    output_["domain"]->open_h5file(true);

  } else {
    for (auto domain = S.mesh_begin(); domain != S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);

      Key domain_name = Keys::replace_all(domain->first, ":", "-");
      std::filesystem::path chkp_file =
        std::filesystem::path(file_or_dirname) / (domain_name + ".h5");
      output_[domain->first] = Teuchos::rcp(new HDF5_MPI(mesh->getComm(), chkp_file.string()));
      output_[domain->first]->open_h5file(true);
    }
  }
}


// Constructor for reading
Checkpoint::Checkpoint(const std::string& filename,
                       const Comm_ptr_type& comm,
                       const std::string& domain)
  : IOEvent()
{
  // if provided a directory, use new style
  if (std::filesystem::is_directory(filename)) {
    Key domain_name = Keys::replace_all(domain, ":", "-");
    std::filesystem::path chkp_file = std::filesystem::path(filename) / (domain_name + ".h5");
    single_file_ = false;
    output_[domain] = Teuchos::rcp(new HDF5_MPI(comm, chkp_file.string()));
    output_[domain]->open_h5file(true);

  } else if (std::filesystem::is_regular_file(filename)) {
    single_file_ = true;
    output_["domain"] = Teuchos::rcp(new HDF5_MPI(comm, filename));
    output_["domain"]->open_h5file(true);
  } else {
    Errors::Message message;
    message << "Checkpoint::Read: location \"" << filename << "\" does not exist.";
    Exceptions::amanzi_throw(message);
  }
}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void
Checkpoint::ReadParameters_()
{
  filebasename_ = plist_.get<std::string>("file name base", "checkpoint");
  filenamedigits_ = plist_.get<int>("file name digits", 5);
  single_file_ = plist_.get<bool>("single file checkpoint", true);
};


void
Checkpoint::CreateFile(const int cycle)
{
  // create the restart directory
  std::stringstream oss;
  oss.flush();
  oss << filebasename_;
  oss.fill('0');
  oss.width(filenamedigits_);
  oss << std::right << cycle;

  if (single_file_) {
    output_["domain"]->createDataFile(oss.str());
    output_["domain"]->open_h5file();

  } else {
    std::filesystem::create_directory(oss.str());
    for (const auto& file_out : output_) {
      Key filename = Keys::replace_all(file_out.first, ":", "-");
      std::filesystem::path chkp_file = std::filesystem::path(oss.str()) / filename;
      file_out.second->createDataFile(chkp_file.string());
      file_out.second->open_h5file();
    }
  }
}


void
Checkpoint::CreateFinalFile(const int cycle)
{
  std::stringstream oss;
  oss.flush();
  oss << filebasename_;
  oss.fill('0');
  oss.width(filenamedigits_);
  oss << std::right << cycle;

  if (single_file_) {
    std::string ch_file = oss.str() + ".h5";
    std::string ch_final = filebasename_ + "_final.h5";
    if (std::filesystem::is_regular_file(ch_final.data())) std::filesystem::remove(ch_final.data());
    std::filesystem::create_hard_link(ch_file.data(), ch_final.data());
  } else {
    std::string ch_dir = oss.str();
    std::string ch_final = filebasename_ + "_final";
    if (std::filesystem::is_directory(ch_final.data())) std::filesystem::remove(ch_final.data());
    std::filesystem::create_symlink(ch_dir.data(), ch_final.data());
  }
}


//
// These are just plain old data in the state now -- can we not read these
// here?  They will get picked up in the data loop of reads.
//
void
Checkpoint::ReadAttributes(State& S)
{}


void
Checkpoint::Finalize()
{
  for (const auto& file_out : output_) { file_out.second->close_h5file(); }
}


void
Checkpoint::Write(const State& S, WriteType write_type, Amanzi::ObservationData* obs_data)
{
  if (!is_disabled()) {
    CreateFile(S.get_cycle());

    // write num procs, as checkpoint files are specific to this
    Write("mpi_num_procs", S.GetMesh()->getComm()->NumProc());

    // create hard link to the final file
    if ((write_type == WriteType::FINAL) && S.GetMesh()->getComm()->MyPID() == 0)
      CreateFinalFile(S.get_cycle());

    for (auto it = S.data_begin(); it != S.data_end(); ++it) {
      it->second->WriteCheckpoint(*this, write_type == WriteType::POST_MORTEM);
    }

    WriteObservations(obs_data);
    Finalize();
  }
}

void
Checkpoint::WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names) const
{
  if (names.size() < vec.NumVectors()) {
    Errors::Message msg("Amanzi::Checkpoint::write_vector... not enough names were specified for "
                        "the the components of the multi vector");
    Exceptions::amanzi_throw(msg);
  }

  if (single_file_) {
    const auto& output = output_.at("domain");
    for (int i = 0; i != vec.NumVectors(); ++i) { output->writeCellDataReal(*vec(i), names[i]); }
  } else {
    // double check that the comms are consistent and that the user is
    // following naming best practices
    Key domain_name = Keys::getDomain(names[0]);
    if (!output_.count(domain_name)) domain_name = "domain";
    const auto& output = output_.at(domain_name);

    if (!sameComm(*output->Comm(), vec.Comm())) {
      Errors::Message msg;
      msg << "Checkpoint::WriteVector : vector \"" << names[0]
          << "\" communicator does not match that of domain \"" << domain_name << "\"";
      Exceptions::amanzi_throw(msg);
    }

    for (int i = 0; i != vec.NumVectors(); ++i) { output->writeCellDataReal(*vec(i), names[i]); }
  }
}


template <>
void
Checkpoint::Write(const std::string& name, const Epetra_Vector& vec) const
{
  if (single_file_) {
    const auto& output = output_.at("domain");
    output->writeCellDataReal(vec, name);
  } else {
    // double check that the comms are consistent and that the user is
    // following naming best practices
    Key domain_name = Keys::getDomain(name);
    if (!output_.count(domain_name)) domain_name = "domain";
    const auto& output = output_.at(domain_name);

    if (!sameComm(*output->Comm(), vec.Comm())) {
      Errors::Message msg;
      msg << "Checkpoint::WriteVector : vector \"" << name
          << "\" communicator does not match that of domain \"" << domain_name << "\"";
      Exceptions::amanzi_throw(msg);
    }

    output->writeCellDataReal(vec, name);
  }
}


// -----------------------------------------------------------------------------
// Write observations
// -----------------------------------------------------------------------------
void
Checkpoint::WriteObservations(ObservationData* obs_data)
{
  if (obs_data == nullptr) return;

  std::vector<std::string> labels = obs_data->observationLabels();
  int nlabels = labels.size();

  int ndata(0);
  double* tmp_data(NULL);

  auto& output = output_["domain"];

  if (nlabels > 0) {
    // save names of observations and their number
    int* nobs;
    char** tmp_labels;

    nobs = (int*)malloc(nlabels * sizeof(int));
    tmp_labels = (char**)malloc(nlabels * sizeof(char*));

    for (int i = 0; i < nlabels; i++) {
      tmp_labels[i] = (char*)malloc((labels[i].size() + 1) * sizeof(char));
      strcpy(tmp_labels[i], labels[i].c_str());
      nobs[i] = (*obs_data)[labels[i]].size();
      ndata += nobs[i];
    }

    output->writeDataString(tmp_labels, nlabels, "obs_names");
    output->writeAttrInt(nobs, nlabels, "obs_numbers");

    // save observation values
    tmp_data = (double*)malloc(2 * ndata * sizeof(double));
    int m(0);
    for (int i = 0; i < nlabels; ++i) {
      std::vector<ObservationData::DataQuadruple> tmp = (*obs_data)[labels[i]];
      for (int k = 0; k < tmp.size(); ++k) {
        tmp_data[m++] = tmp[k].time;
        tmp_data[m++] = tmp[k].value;
      }
    }

    // release memory
    for (int i = 0; i < nlabels; i++) free(tmp_labels[i]);
    free(tmp_labels);
    free(nobs);
  }

  // write data only on the root. This requires global data size.
  int ndata_tmp(ndata);
  output->Comm()->MaxAll(&ndata_tmp, &ndata, 1);

  if (ndata > 0) {
    int ndata_glb = 2 * ndata;
    ndata = (output->Comm()->MyPID() == 0) ? ndata_glb : 0;
    output->writeDatasetReal(tmp_data, ndata, ndata_glb, "obs_values");
    if (tmp_data != NULL) free(tmp_data);
  }
}


} // namespace Amanzi
