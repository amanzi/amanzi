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
#include <iostream>
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "OutputXDMF.hh"
#if ENABLE_Silo
#include "OutputSilo.hh"
#endif

#include "Visualization.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Visualization::Visualization (Teuchos::ParameterList& plist)
  : IOEvent(plist),
    time_unit_written_(false)
{
  ReadParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Visualization  ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_,this);

  my_units_ = plist_.get<std::string>("time units", "y");
  ValidUnitOrThrow_(my_units_);
}


// -----------------------------------------------------------------------------
// Constructor for a disabled Vis.
// -----------------------------------------------------------------------------
Visualization::Visualization ()
  : IOEvent(),
    my_units_("y"),
    time_unit_written_(false)
{}


void Visualization::set_name(const std::string& name) {
  name_ = name;
  domains_.push_back(name);
}

bool Visualization::WritesDomain(const std::string& name) const {
  if (std::find(domains_.begin(), domains_.end(), name) != domains_.end())
    return true;
  if (Keys::isDomainSet(name) && WritesDomain(Keys::getDomainSetName(name)))
    return true;
  return false;
}

void Visualization::AddDomain(const std::string& name) {
  domains_.push_back(name);
}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Visualization::ReadParameters_() {
  Teuchos::Array<std::string> no_regions(0);
  Teuchos::ParameterList& tmp = plist_.sublist("write regions");

  regions_.clear();
  for (Teuchos::ParameterList::ConstIterator it = tmp.begin(); it != tmp.end(); ++it) {
    regions_[it->first] = tmp.get<Teuchos::Array<std::string> >(it->first, no_regions);
  }
  write_partition_ = plist_.get<bool>("write partitions", false);

  dynamic_mesh_ = plist_.get<bool>("dynamic mesh", false);
  write_mesh_exo_ = plist_.get<bool>("write mesh to Exodus II", false);

  if (plist_.isParameter("aliased domains")) {
    auto aliases = plist_.get<Teuchos::Array<std::string>>("aliased domains");
    for (const auto& alias : aliases) {
      domains_.push_back(alias);
    }
  }
}


// -----------------------------------------------------------------------------
// Write a multivector
// -----------------------------------------------------------------------------
void Visualization::WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names) const {
  visualization_output_->WriteMultiVector(vec, names, AmanziMesh::CELL);
}


// -----------------------------------------------------------------------------
// Write a vector
// -----------------------------------------------------------------------------
void Visualization::WriteVector(const Epetra_Vector& vec, const std::string& name) const {
  visualization_output_->WriteVector(vec, name, AmanziMesh::CELL);
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
          mesh_->get_set_entities(*reg_it, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &ids);

          for (auto jt = ids.begin(); jt != ids.end(); ++jt) {
            reg[0][*jt] = reg_index;
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
void Visualization::CreateFiles(bool include_io_set) {
  AMANZI_ASSERT(mesh_ != Teuchos::null);

  std::string file_format = plist_.get<std::string>("file format", "XDMF");

  if (file_format == "XDMF" || file_format == "xdmf") {
    visualization_output_ = Teuchos::rcp(new OutputXDMF(plist_, mesh_, true, dynamic_mesh_, include_io_set));
#if ENABLE_Silo
  } else if (file_format == "Silo" || file_format == "SILO" || file_format == "silo") {
    visualization_output_ = Teuchos::rcp(new OutputSilo(plist_, mesh_, true, dynamic_mesh_));
#endif
  } else {
    Errors::Message msg("Visualization: Unknown file format: \""+file_format+"\"");
    Exceptions::amanzi_throw(msg);
  }

  if (write_mesh_exo_ && !dynamic_mesh_) {
    std::stringstream mesh_fname;
    mesh_fname << plist_.get<std::string>("file name base") << "_mesh.exo";
    mesh_->write_to_exodus_file(mesh_fname.str());
  }
}


void Visualization::CreateTimestep(double time, int cycle, const std::string& tag) {
  bool success = false;
  time = units_.ConvertTime(time, "s", my_units_, success);
  AMANZI_ASSERT(success);

  visualization_output_->InitializeCycle(time, cycle, tag);
  if (!time_unit_written_) {
    visualization_output_->WriteAttribute(my_units_, "time unit");
    time_unit_written_ = true;
  }

  if (write_mesh_exo_ && dynamic_mesh_) {
    std::stringstream mesh_fname;
    mesh_fname << plist_.get<std::string>("file name base") << "_mesh_" << cycle << ".exo";
    mesh_->write_to_exodus_file(mesh_fname.str());
  }
}


void Visualization::FinalizeTimestep() const {
  visualization_output_->FinalizeCycle();
}

} // namespace Amanzi
