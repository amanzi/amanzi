/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Visualization of data.

------------------------------------------------------------------------- */

#include "visualization.hh"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Visualization::Visualization (Teuchos::ParameterList& plist, Epetra_MpiComm* comm) :
  IOEvent(plist, comm), dynamic_mesh_(false) {
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Visualization  ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_,this);

  // Set up the HDF5
  visualization_output_ = Teuchos::rcp(new Amanzi::HDF5_MPI(*comm_));
  visualization_output_->setTrackXdmf(true);
  visualization_output_->setDynMesh(dynamic_mesh_);
}


// -----------------------------------------------------------------------------
// Constructor for a disabled Vis.
// -----------------------------------------------------------------------------
Visualization::Visualization () : IOEvent() {}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Visualization::ReadParameters_() {
  filebasename_ = plist_.get<std::string>("file name base","amanzi_vis");
  dynamic_mesh_ = plist_.get<bool>("dynamic mesh",false);
}


// -----------------------------------------------------------------------------
// Write a multivector
// -----------------------------------------------------------------------------
void Visualization::WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const {
  if (names.size() < vec.NumVectors()) {
    Errors::Message m("Amanzi::Visualization::write_vector... not enough names were specified for the the components of the multi vector");
    Exceptions::amanzi_throw(m);
  }
  for (int i=0; i!=vec.NumVectors(); ++i) {
    visualization_output_->writeCellDataReal(*vec(i), names[i]);
  }
}


// -----------------------------------------------------------------------------
// Write a vector
// -----------------------------------------------------------------------------
void Visualization::WriteVector(const Epetra_Vector& vec, const std::string name ) const {
  visualization_output_->writeCellDataReal(vec ,name);
}



// -----------------------------------------------------------------------------
// Write the mesh
// -----------------------------------------------------------------------------
void Visualization::WriteMesh(const double time, const int iteration) const {
  visualization_output_->writeMesh(time, iteration);
}


void Visualization::CreateFiles() {

  if (!is_disabled()) {
    // create file name for the mesh
    std::stringstream meshfilename;
    meshfilename.flush();
    meshfilename << filebasename_;
    meshfilename << "_mesh";
    // create file name for the data
    std::stringstream datafilename;
    datafilename.flush();
    datafilename << filebasename_;
    datafilename << "_data";
    // create the files
    visualization_output_->createMeshFile(mesh_, meshfilename.str());
    visualization_output_->createDataFile(datafilename.str());
  }
}

void Visualization::CreateTimestep(const double& time, const int& cycle) {
  visualization_output_->createTimestep(time,cycle);
}

void Visualization::FinalizeTimestep() const {
  visualization_output_->endTimestep();
}

} // namespace
