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
#include <iostream>
#include <iomanip>

#include "boost/filesystem.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Output.hh"
#include "OutputFactory.hh"
#include "Input.hh"
#include "InputFactory.hh"
#include "State.hh"
#include "Checkpoint.hh"

namespace Amanzi {

Checkpoint::Checkpoint(bool single_file) : single_file_(single_file)
{
  Teuchos::ParameterList plist("domain");
  output_["domain"] = OutputFactory::createForCheckpoint(plist, getDefaultComm());
}

Checkpoint::Checkpoint(Teuchos::ParameterList& plist, const State& S)
  : IOEvent(plist),
    single_file_(plist.get<bool>("single file", true))
{
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Checkpoint     ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_, this);

  // get the file name base
  filenamebase_ = plist.get<std::string>("file name base", "checkpoint");

  // Set up the HDF5 objects
  if (single_file_) {
    // Single_File style is one checkpoint file on "domain"'s comm
    Comm_ptr_type comm;
    if (S.HasMesh("domain")) {
      comm = S.GetMesh()->getComm();
    } else {
      comm = Amanzi::getDefaultComm();
    }
    plist.setName("domain");
    output_["domain"] = OutputFactory::createForCheckpoint(plist, comm);

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
    plist.set<bool>("single file", false);
    for (auto domain = S.mesh_begin(); domain != S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);
      plist.setName(domain->first);
      output_[domain->first] = OutputFactory::createForCheckpoint(plist, mesh->getComm());
    }
  }
}


// Constructor for reading
Checkpoint::Checkpoint(const std::string& file_or_dirname, const State& S)
  : IOEvent()
{
  // if provided a directory, use new style
  if (boost::filesystem::is_directory(file_or_dirname)) {
    single_file_ = false;
  } else if (boost::filesystem::is_regular_file(file_or_dirname)) {
    single_file_ = true;
  } else {
    Errors::Message message;
    message << "Checkpoint::Read: location \"" << file_or_dirname << "\" does not exist.";
    Exceptions::amanzi_throw(message);
  }

  // create the readers
  if (single_file_) {
    Teuchos::ParameterList plist("checkpoint");
    plist.set<std::string>("file name", file_or_dirname);
    auto comm = S.GetMesh()->getComm();
    input_["domain"] = InputFactory::createForCheckpoint(plist, comm);

  } else {
    for (auto domain = S.mesh_begin(); domain != S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);

      Key domain_name = Keys::replace_all(domain->first, ":", "-");
      boost::filesystem::path chkp_file =
        boost::filesystem::path(file_or_dirname) / (domain_name + ".h5");
      Teuchos::ParameterList plist(domain_name);
      plist.set<std::string>("file name", chkp_file.string());
      input_[domain->first] = InputFactory::createForCheckpoint(plist, mesh->getComm());
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
  if (boost::filesystem::is_directory(filename)) {
    single_file_ = false;

    Key domain_name = Keys::replace_all(domain, ":", "-");
    boost::filesystem::path chkp_file = boost::filesystem::path(filename) / (domain_name + ".h5");

    Teuchos::ParameterList plist(domain_name);
    plist.set<std::string>("file name", chkp_file.string());
    plist.set<std::string>("file format", "HDF5");
    filenamebase_ = filename;
    input_[domain] = InputFactory::createForCheckpoint(plist, comm);

  } else if (boost::filesystem::is_regular_file(filename)) {
    Teuchos::ParameterList plist("input checkpoint");
    plist.set<std::string>("file name", filename);
    plist.set<std::string>("file format", "HDF5");
    filenamebase_ = filename;
    single_file_ = true;
    input_["domain"] = InputFactory::createForCheckpoint(plist, comm);

  } else {
    Errors::Message message;
    message << "Checkpoint::Read: location \"" << filename << "\" does not exist.";
    Exceptions::amanzi_throw(message);
  }

  // check the comm_size for consistency
  int file_num_procs(-1);
  read(Teuchos::ParameterList("mpi_num_procs"), file_num_procs);
  if (comm->getSize() != file_num_procs) {
    Errors::Message msg;
    msg << "Requested checkpoint file " << filename << " was created on " << file_num_procs
        << " processes, making it incompatible with this run on " << comm->getSize()
        << " processes.";
    Exceptions::amanzi_throw(msg);
  }
}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void
Checkpoint::readParameters_()
{
  single_file_ = plist_.get<bool>("single file checkpoint", true);
};


void
Checkpoint::createFile_(const int cycle)
{
  if (single_file_) {
    output_["domain"]->createTimestep(-1.0, cycle);

  } else {
    std::string dirname = output_["domain"]->getFilename(cycle);
    boost::filesystem::create_directory(dirname);
    for (const auto& file_out : output_) {
      file_out.second->createTimestep(-1.0, cycle);
    }
  }
}


void
Checkpoint::createFinalFile_(int cycle)
{
  if (single_file_) {
    std::string ch_file = output_["domain"]->getFilename(cycle);
    std::string ch_final = filenamebase_ + "_final.h5";
    if (boost::filesystem::is_regular_file(ch_final.data()))
      boost::filesystem::remove(ch_final.data());
    boost::filesystem::create_hard_link(ch_file.data(), ch_final.data());
  } else {
    std::string ch_file = output_["domain"]->getFilename(cycle);
    std::string ch_final = filenamebase_ + "_final";
    if (boost::filesystem::is_directory(ch_final.data()))
      boost::filesystem::remove(ch_final.data());
    boost::filesystem::create_symlink(ch_final.data(), ch_final.data());
  }
}


void
Checkpoint::finalizeFile_()
{
  for (const auto& file_out : output_) { file_out.second->finalizeTimestep(); }
}


void
Checkpoint::write(const State& S, WriteType write_type, Amanzi::ObservationData* obs_data)
{
  if (!is_disabled()) {
    createFile_(S.get_cycle());

    // write num procs, as checkpoint files are specific to this
    write(Teuchos::ParameterList("mpi_num_procs"), S.GetMesh()->getComm()->getSize());

    // write the data
    for (auto it = S.data_begin(); it != S.data_end(); ++it) {
      it->second->WriteCheckpoint(*this, write_type == WriteType::POST_MORTEM);
    }

    writeObservations_(obs_data);
    finalizeFile_();

    // create hard link to the final file
    if ((write_type == WriteType::FINAL) && S.GetMesh()->getComm()->getRank() == 0)
      createFinalFile_(S.get_cycle());
  }
}


// -----------------------------------------------------------------------------
// Write observations
// -----------------------------------------------------------------------------
void
Checkpoint::writeObservations_(ObservationData* obs_data)
{
  // if (obs_data == nullptr) return;

  // Teuchos::Array<std::string> labels = obs_data->observationLabels();
  // int nlabels = labels.size();

  // int ndata(0);
  // double* tmp_data(NULL);

  // auto& output = output_["domain"];

  // if (nlabels > 0) {
  //   // save names of observations and their number
  //   int* nobs;
  //   char** tmp_labels;

  //   nobs = (int*)malloc(nlabels * sizeof(int));
  //   tmp_labels = (char**)malloc(nlabels * sizeof(char*));

  //   for (int i = 0; i < nlabels; i++) {
  //     tmp_labels[i] = (char*)malloc((labels[i].size() + 1) * sizeof(char));
  //     strcpy(tmp_labels[i], labels[i].c_str());
  //     nobs[i] = (*obs_data)[labels[i]].size();
  //     ndata += nobs[i];
  //   }

  //   Teuchos::ParameterList attrs1("obs_names");
  //   output->write(tmp_labels, nlabels, "obs_names");
  //   output->writeAttrInt(nobs, nlabels, "obs_numbers");

  //   // save observation values
  //   tmp_data = (double*)malloc(2 * ndata * sizeof(double));
  //   int m(0);
  //   for (int i = 0; i < nlabels; ++i) {
  //     std::vector<ObservationData::DataQuadruple> tmp = (*obs_data)[labels[i]];
  //     for (int k = 0; k < tmp.size(); ++k) {
  //       tmp_data[m++] = tmp[k].time;
  //       tmp_data[m++] = tmp[k].value;
  //     }
  //   }

  //   // release memory
  //   for (int i = 0; i < nlabels; i++) free(tmp_labels[i]);
  //   free(tmp_labels);
  //   free(nobs);
  // }

  // // write data only on the root. This requires global data size.
  // int ndata_tmp(ndata);
  // output->Comm()->MaxAll(&ndata_tmp, &ndata, 1);

  // if (ndata > 0) {
  //   int ndata_glb = 2 * ndata;
  //   ndata = (output->Comm()->MyPID() == 0) ? ndata_glb : 0;
  //   output->writeDatasetReal(tmp_data, ndata, ndata_glb, "obs_values");
  //   if (tmp_data != NULL) free(tmp_data);
  // }
}


void
Checkpoint::read(State& S)
{
  // Load the number of processes and ensure they are the same.
  int file_num_procs(-1);
  read(Teuchos::ParameterList("mpi_num_procs"), file_num_procs);
  int my_num_procs = S.GetMesh("domain")->getComm()->getSize();
  if (my_num_procs != file_num_procs) {
    Errors::Message msg;
    msg << "Requested checkpoint file " << filenamebase_ << " was created on " << file_num_procs
        << " processes, making it incompatible with this run on " << my_num_procs
        << " processes.";
    throw(msg);
  }

  // create other files if not single-file
  if (!single_file_) {
    for (auto domain = S.mesh_begin(); domain != S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);

      Key domain_name = Keys::cleanName(Keys::replace_all(domain->first, ":", "-"));
      boost::filesystem::path chkp_file = boost::filesystem::path(filenamebase_) / (domain_name + ".h5");

      Teuchos::ParameterList plist(domain_name);
      plist.set<std::string>("file name", chkp_file.string());
      plist.set<std::string>("file format", "HDF5");

      input_[domain->first] = InputFactory::createForCheckpoint(plist, mesh->getComm());
    }
  }

  // load the data
  for (auto data = S.data_begin(); data != S.data_end(); ++data) {
    data->second->ReadCheckpoint(*this);
  }
}


} // namespace Amanzi
