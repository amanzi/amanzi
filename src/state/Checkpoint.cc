/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Checkpointing for state.
------------------------------------------------------------------------- */
#include <iostream>
#include <iomanip>

#include "mpi.h"
#include "boost/filesystem.hpp"

#include "Checkpoint.hh"
#include "State.hh"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Amanzi {

Checkpoint::Checkpoint(Teuchos::ParameterList& plist,
                       const State& S) :
  IOEvent(plist),
  old_(true)
{
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Checkpoint     ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_,this);

  // Set up the HDF5 objects
  if (old_) {
    // Old style is one checkpoint file on MPI_COMM_WORLD
    auto comm = Amanzi::getDefaultComm();
    output_["domain"] = Teuchos::rcp(new HDF5_MPI(comm));
    output_["domain"]->setTrackXdmf(false);

    // This requires that all meshes and data are defined on MPI_COMM_WORLD
    // Check to confirm!
    for (auto entry=S.mesh_begin(); entry!=S.mesh_end(); ++entry) {
      const std::string& domain = entry->first;
      const auto& mesh = S.GetMesh(domain);
      if (!sameComm(*comm, *mesh->get_comm())) {
        std::stringstream msg;
        msg << "Checkpointing: cannot use single file checkpointing when not all meshes are on MPI_COMM_WORLD (mesh \""
            << domain << "\" not on MPI_COMM_WORLD).  Using multi-file checkpointing.  (Hide this warning by setting \"single file checkpoint\" to \"false\" in the \"checkpoint\" list.";
        std::cerr << "WARNING: " << msg.str() << std::endl;
        old_ = false;
        break;
      }
    }
  }

  // NOTE: do not make this an 'else' clause!
  if (!old_) {
    for (auto domain=S.mesh_begin(); domain!=S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);
      output_[domain->first] = Teuchos::rcp(new HDF5_MPI(mesh->get_comm()));
      output_[domain->first]->setTrackXdmf(false);
    }
  }
}


// this constructor makes an object for reading
Checkpoint::Checkpoint(bool old) :
  IOEvent(),
  old_(old)
{}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Checkpoint::ReadParameters_() {
  filebasename_ = plist_.get<std::string>("file name base","checkpoint");
  filenamedigits_ = plist_.get<int>("file name digits", 5);
  old_ = plist_.get<bool>("single file checkpoint", true);
};


void Checkpoint::CreateFile(const int cycle) {
  // create the restart directory
  std::stringstream oss;
  oss.flush();
  oss << filebasename_;
  oss.fill('0');
  oss.width(filenamedigits_);
  oss << std::right << cycle;

  if (old_) {
    output_["domain"]->createDataFile(oss.str());
    output_["domain"]->open_h5file();

  } else {
    boost::filesystem::create_directory(oss.str());
    for (const auto& file_out : output_) {
      boost::filesystem::path chkp_file = boost::filesystem::path(oss.str()) / file_out.first;
      file_out.second->createDataFile(chkp_file.string());
      file_out.second->open_h5file();
    }
  }
}


void Checkpoint::CreateFinalFile(const int cycle) {
  std::stringstream oss;
  oss.flush();
  oss << filebasename_;
  oss.fill('0');
  oss.width(filenamedigits_);
  oss << std::right << cycle;

  if (old_) {
    std::string ch_file = oss.str() + ".h5";
    std::string ch_final = filebasename_ + "_final.h5";
    if (boost::filesystem::is_regular_file(ch_final.data()))
      boost::filesystem::remove(ch_final.data());
    boost::filesystem::create_hard_link(ch_file.data(), ch_final.data());
  } else {

    std::string ch_dir = oss.str();
    std::string ch_final = filebasename_ + "_final";
    if (boost::filesystem::is_directory(ch_final.data()))
      boost::filesystem::remove(ch_final.data());
    boost::filesystem::create_symlink(ch_dir.data(), ch_final.data());
  }
}


void Checkpoint::Finalize() {
  for (const auto& file_out : output_) {
    file_out.second->close_h5file();
  }
}

void Checkpoint::Write(const State& S,
                       double dt,
                       bool final,
                       Amanzi::ObservationData* obs_data)
{
  if (!is_disabled()) {
    CreateFile(S.cycle());

    // create hard link to the final file
    if (final && S.GetMesh()->get_comm()->MyPID() == 0)
      CreateFinalFile(S.cycle());

    for (State::field_iterator field=S.field_begin(); field!=S.field_end(); ++field) {
      field->second->WriteCheckpoint(*this);
    }
    WriteAttributes(S.GetMesh("domain")->get_comm()->NumProc(),
                    S.time(), dt, S.cycle(), S.position());
    WriteObservations(obs_data);
    Finalize();
  }
}

double Checkpoint::Read(State& S, const std::string& file_or_dirname)
{
  // create the readers
  auto comm = S.GetMesh()->get_comm();
  if (old_) {
    output_["domain"] = Teuchos::rcp(new HDF5_MPI(comm, file_or_dirname));
    output_["domain"]->open_h5file();

  } else {
    for (auto domain=S.mesh_begin(); domain!=S.mesh_end(); ++domain) {
      const auto& mesh = S.GetMesh(domain->first);

      boost::filesystem::path chkp_file = boost::filesystem::path(file_or_dirname) / (domain->first+".h5");
      output_[domain->first] = Teuchos::rcp(new HDF5_MPI(mesh->get_comm(), chkp_file.string()));
      output_[domain->first]->open_h5file();
    }
  }

  // load the attributes
  double time(0.);
  output_["domain"]->readAttrReal(time, "time");
  S.set_time(time);

  double dt(0.);
  output_["domain"]->readAttrReal(dt, "dt");

  int cycle(0);
  output_["domain"]->readAttrInt(cycle, "cycle");
  S.set_cycle(cycle);

  int pos(0);
  output_["domain"]->readAttrInt(pos, "position");
  S.set_position(pos);

  // load the number of processes and ensure they are the same -- otherwise
  // the below just gives crap.
  int rank(-1);
  output_["domain"]->readAttrInt(rank, "mpi_comm_world_rank");
  if (comm->NumProc() != rank) {
    std::stringstream messagestream;
    messagestream << "Requested checkpoint file " << file_or_dirname << " was created on "
                  << rank << " processes, making it incompatible with this run on "
                  << comm->NumProc() << " process.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // load the data
  for (State::field_iterator field=S.field_begin(); field!=S.field_end(); ++field) {
    if (field->second->type() == COMPOSITE_VECTOR_FIELD &&
        field->second->io_checkpoint()) {
      Key domain = Keys::getDomain(field->first);
      if (domain.empty()) domain = "domain";

      bool read_complete = field->second->ReadCheckpoint(*output_[domain]);
      if (read_complete) field->second->set_initialized();
    }
  }

  // close files
  for (auto& file_out : output_) {
    file_out.second->close_h5file();
  }
  return dt;
}



void Checkpoint::WriteVector(const Epetra_MultiVector& vec,
                             const std::vector<std::string>& names ) const {
  if (names.size() < vec.NumVectors()) {
    Errors::Message msg("Amanzi::Checkpoint::write_vector... not enough names were specified for the the components of the multi vector");
    Exceptions::amanzi_throw(msg);
  }

  if (old_) {
    const auto& output = output_.at("domain");
    for (int i=0; i!=vec.NumVectors(); ++i) {
      output->writeCellDataReal(*vec(i), names[i]);
    }
  } else {
    // double check that the comms are consistent and that the user is
    // following naming best practices
    Key domain_name = Keys::getDomain(names[0]);
    if (domain_name.empty()) domain_name = "domain";
    const auto& output = output_.at(domain_name);

    if (!sameComm(*output->Comm(), vec.Comm())) {
      Errors::Message msg;
      msg << "Checkpoint::WriteVector : vector \"" << names[0] << "\" communicator does not match that of domain \"" << domain_name << "\"";
      Exceptions::amanzi_throw(msg);
    }

    for (int i=0; i!=vec.NumVectors(); ++i) {
      output->writeCellDataReal(*vec(i), names[i]);
    }
  }
}


// -----------------------------------------------------------------------------
// Write simple attributes.
// -----------------------------------------------------------------------------
void Checkpoint::WriteAttributes(int comm_size, double time, double dt, int cycle, int position) const {
  const auto& output = output_.at("domain");
  output->writeAttrReal(time, "time");
  output->writeAttrReal(dt, "dt");
  output->writeAttrInt(cycle, "cycle");
  output->writeAttrInt(position, "position");
  output->writeAttrInt(comm_size, "mpi_comm_world_rank");
};


void Checkpoint::WriteAttributes(int comm_size, double time, double dt, int cycle) const {
  const auto& output = output_.at("domain");
  output->writeAttrReal(time, "time");
  output->writeAttrReal(dt, "dt");
  output->writeAttrInt(cycle, "cycle");
  output->writeAttrInt(comm_size, "mpi_comm_world_rank");
};


void Checkpoint::WriteAttributes(int comm_size, double time, int cycle) const {
  const auto& output = output_.at("domain");
  output->writeAttrReal(time, "time");
  output->writeAttrInt(cycle, "cycle");
  output->writeAttrInt(comm_size, "mpi_comm_world_rank");
};


// -----------------------------------------------------------------------------
// Write observations
// -----------------------------------------------------------------------------
void Checkpoint::WriteObservations(ObservationData* obs_data)
{
  if (obs_data == nullptr) return;

  std::vector<std::string> labels = obs_data->observationLabels();
  int nlabels = labels.size();

  int ndata(0);
  double* tmp_data(NULL);

  auto& output = output_["domain"];

  if (nlabels > 0) {
    // save names of observations and their number
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

  // write data only on the root
  if (ndata > 0) {
    int ndata_glb = 2 * ndata;
    ndata = (output->Comm()->MyPID() == 0) ? ndata_glb : 0;
    output->writeDatasetReal(tmp_data, ndata, ndata_glb, "obs_values");
    if (tmp_data != NULL) free(tmp_data);
  }
}

}  // namespace Amanzi
