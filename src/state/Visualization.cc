/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)

  Visualization of data.
*/

// TPLs
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

// State
#include "Visualization.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Visualization::Visualization (Teuchos::ParameterList& plist, Epetra_MpiComm* comm) :
  IOEvent(plist), dynamic_mesh_(false) {
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Visualization  ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_,this);

  // Set up the HDF5
  visualization_output_ = Teuchos::rcp(new Amanzi::HDF5_MPI(*comm));
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

  Teuchos::Array<std::string> no_regions(0);
  Teuchos::ParameterList& tmp = plist_.sublist("write regions");
    
  regions_.clear();
  for (Teuchos::ParameterList::ConstIterator it = tmp.begin(); it != tmp.end(); ++it) {
    regions_[it->first] = tmp.get<Teuchos::Array<std::string> >(it->first, no_regions);
  }
  write_partition_ = plist_.get<bool>("write partitions", false);
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


// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
void Visualization::WriteRegions() {
  if (regions_.size() > 0) {
    for (std::map<std::string, Teuchos::Array<std::string> >::const_iterator it = regions_.begin();
         it != regions_.end(); ++it) {
      // first make an Epetra_Vector to hold the region information
      Epetra_MultiVector reg(mesh_->cell_map(false), 1, true);

      // loop over the regions and initialize the reg array
      double reg_index = 1.0;
      for (Teuchos::Array<std::string>::const_iterator reg_it = (it->second).begin();
         reg_it != (it->second).end(); ++reg_it, reg_index += 1.0) {
        // only do something if the user provided a valid region name
        // for a region that consists of cells
        if (mesh_->valid_set_name(*reg_it, AmanziMesh::CELL)) {
          AmanziMesh::Entity_ID_List ids;
          mesh_->get_set_entities(*reg_it, AmanziMesh::CELL, AmanziMesh::OWNED, &ids);

          for (AmanziMesh::Entity_ID_List::const_iterator it = ids.begin(); it != ids.end(); ++it) {
            reg[0][*it] = reg_index;
          }
        }
      }
      std::vector<std::string> name; 
      name.push_back(it->first);
      WriteVector(reg, name);
    }
  }
}


// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
void Visualization::WritePartition() {
  if (write_partition_) {
    // first make an Epetra_Vector to hold the partitions information
    Epetra_MultiVector reg(mesh_->cell_map(false),1,false);
    // loop over the regions and initialize the reg array
    double part_index = static_cast<double>(mesh_->get_comm()->MyPID());
    reg.PutScalar(part_index);
  
    std::vector<std::string> name; 
    name.push_back("partition");
    WriteVector(reg,name);
  }
}


// -----------------------------------------------------------------------------
// Writing to files
// -----------------------------------------------------------------------------
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
  visualization_output_->open_h5file();
}


void Visualization::FinalizeTimestep() const {
  visualization_output_->close_h5file();
  visualization_output_->endTimestep();
}

} // namespace Amanzi
