/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
      Markus Berndt
*/

//! Checkpoint writes ALL data, but no meshes to files.

/*
  Reads/writes to/from file using a generic Input/Output object.
*/

#include <iomanip>
#include <iostream>

#include "boost/filesystem.hpp"

#include "Checkpoint.hh"
#include "UniqueHelpers.hh"
#include "InputFactory.hh"
#include "OutputFactory.hh"

namespace Amanzi {

// standard constructor
Checkpoint::Checkpoint(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                       const Comm_ptr_type& comm,
                       bool read)
  : IOEvent(plist),
    filename_base_(
      plist->template get<std::string>("file name base", "checkpoint")),
    comm_(comm)

{
  if (read)
    input_ = std::move(InputFactory::CreateForCheckpoint(*plist, comm));
  else
    output_ = std::move(OutputFactory::CreateForCheckpoint(*plist, comm));
}

void
Checkpoint::CreateFile(double time, int cycle)
{
  output_->CreateFile(time, cycle);

  // write comm size
  Teuchos::ParameterList plist("mpi_num_procs");
  plist.set("units", "-");
  Write(plist, comm_->getSize());
};

void
Checkpoint::FinalizeFile(bool final)
{
  output_->FinalizeFile();
  if (final) {
    std::string filename = output_->Filename();
    std::string filename_extension = boost::filesystem::extension(filename);
    std::string filename_final =
      filename_base_ + "_finale" + filename_extension;
    boost::filesystem::remove(filename_extension);
    boost::filesystem::create_hard_link(filename.data(), filename_final.data());
  }
}

} // namespace Amanzi
