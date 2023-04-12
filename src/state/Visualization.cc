/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  State

  Visualization of data.

*/

// TPLs
#include <iostream>
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "OutputFactory.hh"
#include "State.hh"
#include "Visualization.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Visualization::Visualization(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             bool include_io_set)
  : IOEvent(plist),
    time_unit_written_(false),
    count_(0),
    mesh_(mesh)
{
  readParameters_();

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Visualization  ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist_, this);

  my_units_ = plist_.get<std::string>("time units", "y");
  ValidUnitOrThrow_(my_units_);

  // if we are not including io set, there is only one set of files -- create
  // them now and this is RAII
  if (!include_io_set && mesh_ != Teuchos::null) createFiles(include_io_set);
}


// -----------------------------------------------------------------------------
// Constructor for a disabled Vis.
// -----------------------------------------------------------------------------
Visualization::Visualization()
  : IOEvent(),
    my_units_("y"),
    time_unit_written_(false),
    count_(0) {}


bool
Visualization::writesDomain(const std::string& name) const
{
  if (std::find(domains_.begin(), domains_.end(), name) != domains_.end()) return true;
  if (Keys::isDomainSet(name) && writesDomain(Keys::getDomainSetName(name))) return true;
  return false;
}

void
Visualization::addDomain(const std::string& name)
{
  domains_.push_back(name);
}


// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void
Visualization::readParameters_()
{
  addDomain(Keys::cleanPListName(plist_));

  Teuchos::Array<std::string> no_regions(0);
  Teuchos::ParameterList& tmp = plist_.sublist("write regions");

  regions_.clear();
  for (Teuchos::ParameterList::ConstIterator it = tmp.begin(); it != tmp.end(); ++it) {
    regions_[it->first] = tmp.get<Teuchos::Array<std::string>>(it->first, no_regions);
  }
  write_partition_ = plist_.get<bool>("write partitions", false);

  dynamic_mesh_ = plist_.get<bool>("dynamic mesh", false);
  write_mesh_exo_ = plist_.get<bool>("write mesh to Exodus II", false);

  if (plist_.isParameter("aliased domains")) {
    auto aliases = plist_.get<Teuchos::Array<std::string>>("aliased domains");
    for (const auto& alias : aliases) {
      KeyTriple ds;
      if (Keys::splitDomainSet(alias, ds)) {
        domains_.push_back(std::get<0>(ds));
      } else {
        domains_.push_back(alias);
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
void
Visualization::writeRegions() const
{
  if (regions_.size() > 0) {
    for (std::map<std::string, Teuchos::Array<std::string>>::const_iterator it = regions_.begin();
         it != regions_.end();
         ++it) {
      // first make an Epetra_Vector to hold the region information
      IntVector_type reg(mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1);

      // loop over the regions and initialize the reg array
      {
        auto reg_v = reg.getLocalViewHost(Tpetra::Access::ReadWrite);
        double reg_index = 1.0;
        for (Teuchos::Array<std::string>::const_iterator reg_it = (it->second).begin();
             reg_it != (it->second).end();
             ++reg_it, reg_index += 1.0) {
          // only do something if the user provided a valid region name
          // for a region that consists of cells
          if (mesh_->isValidSetName(*reg_it, AmanziMesh::Entity_kind::CELL)) {
            auto ids = mesh_->getSetEntities(
              *reg_it, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

            for (auto jt : ids) reg_v(jt, 0) = reg_index;
          }
        }
      }

      Teuchos::ParameterList attrs(it->first);
      attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);
      write(attrs, reg);
    }
  }
}


// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
void
Visualization::writePartition() const
{
  if (write_partition_) {
    // first make an Epetra_Vector to hold the partitions information
    IntVector_type reg(mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
    // loop over the regions and initialize the reg array
    int part_index = mesh_->getComm()->getRank();
    reg.putScalar(part_index);

    Teuchos::ParameterList attrs("partition");
    attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);
    write(attrs, reg);
  }
}


// -----------------------------------------------------------------------------
// Writing to files
// -----------------------------------------------------------------------------
void
Visualization::createFiles(bool include_io_set)
{
  AMANZI_ASSERT(mesh_ != Teuchos::null);

  Teuchos::ParameterList output_plist(plist_);
  if (include_io_set)
    output_plist.set("file name base", Keys::cleanPListName(plist_) + "-io-counter" + std::to_string(count_++));
  output_ = OutputFactory::createForVis(output_plist, AmanziMesh::onMemSpace<MemSpace_kind::HOST>(mesh_));

  if (write_mesh_exo_ && !dynamic_mesh_ && mesh_->getMeshFramework()) {
    std::string mesh_fname = output_plist.get<std::string>("file name base") + "_mesh.exo";
    mesh_->getMeshFramework()->writeToExodusFile(mesh_fname);
  }
}


void
Visualization::createTimestep(double time, int cycle)
{
  bool success = false;
  time = units_.ConvertTime(time, "s", my_units_, success);
  AMANZI_ASSERT(success);

  output_->createTimestep(time, cycle);
  if (!time_unit_written_) {
    output_->write(Teuchos::ParameterList("time unit"), my_units_);
    time_unit_written_ = true;
  }

  if (write_mesh_exo_ && dynamic_mesh_) {
    std::stringstream mesh_fname;
    mesh_fname << plist_.get<std::string>("file name base") << "_mesh_" << cycle << ".exo";
    mesh_->getMeshFramework()->writeToExodusFile(mesh_fname.str());
  }
}


void
Visualization::finalizeTimestep()
{
  output_->finalizeTimestep();
}


void
Visualization::write(const State& S)
{
  if (!is_disabled()) {
    // Create the new time step
    createTimestep(S.get_time(), S.get_cycle());

    for (auto r = S.data_begin(); r != S.data_end(); ++r) {
      if (writesDomain(Keys::getDomain(r->first)) ||
          r->first == "time" ||
          r->first == "cycle") {
        // Should we vis all tags or just the default tag?
        // -- write all tags
        // r->WriteVis(vis, nullptr);

        // -- write default tag
        // Tag tag;
        // r->second->WriteVis(vis, &tag);

        // -- write default tag if it exists, else write another tag with the
        // -- same time
        Tag tag;
        if (r->second->HasRecord(tag)) {
          r->second->WriteVis(*this, &tag);
        } else {
          // try to find a record at the same time
          double time = S.get_time();
          for (const auto& time_record : S.GetRecordSet("time")) {
            if (r->second->HasRecord(time_record.first) &&
                S.get_time(time_record.first) == time) {
              r->second->WriteVis(*this, &time_record.first);
              break;
            }
          }
        }
      }
    }
    writeRegions();
    writePartition();
    finalizeTimestep();
  }
}

} // namespace Amanzi
