/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Factory for Output objects.
/*

*/

#include <string>

#include "errors.hh"
#include "Key.hh"

#include "Output.hh"
#include "OutputXDMF.hh"
#include "InputOutputHDF5.hh"
#include "InputOutputVerticallyStructuredNetCDF.hh"
#include "IOFileFormat.hh"
#ifdef HAVE_SILO
#  include "OutputSilo.hh"
#endif

#include "OutputFactory.hh"


namespace Amanzi {
namespace OutputFactory {

// creates a formatter for the directory name of a multi-file checkpoint/vis
Output::FilenameFormatter
createDirectoryFormatter(Teuchos::ParameterList& plist)
{
  std::string filename_base = plist.get<std::string>("file name base");
  int sig_digs = plist.get<int>("file name base digits", 5);

  auto formatter = [=](int count) {
    std::stringstream str;
    str << filename_base;
    str.fill('0');
    str.width(sig_digs);
    str << std::right << count;
    return str.str();
  };
  return formatter;
}


std::unique_ptr<Output>
createForVis(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  // set the filename base based on domain if needed
  std::string domain_name = Keys::cleanName(Keys::cleanPListName(plist));
  std::string filenamebase = plist.get<std::string>("file name base", "amanzi_vis");
  if (!domain_name.empty() && domain_name != "domain")
    filenamebase = filenamebase + "_" + domain_name;
  plist.set<std::string>("file name base", filenamebase);

  Utils::IOFileFormat_kind file_format =
    Utils::createIOFileFormat(plist, Utils::IOFileFormat_kind::XDMF);
  switch (file_format) {
  case Utils::IOFileFormat_kind::XDMF:
    return std::make_unique<OutputXDMF>(plist, mesh);
  case Utils::IOFileFormat_kind::NETCDF:
    return std::make_unique<OutputVerticallyStructuredNetCDF>(plist, mesh->getComm(), mesh);
  case Utils::IOFileFormat_kind::SILO:
#ifdef HAVE_SILO
    return std::make_unique<OutputSilo>(plist, mesh);
#else
    {
      Errors::Message msg(
        "OutputFactory: requested output type Silo was not enabled at configure time.");
      throw(msg);
    }
#endif
  default: {
    Errors::Message msg;
    msg << "OutputFactory::createForVis: invalid \"file format\" \""
        << Utils::to_string(file_format) << "\" for writing vis";
    throw(msg);
  }
  }
}

std::unique_ptr<Output>
createForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm)
{
  // set the filename base based on domain if needed
  std::string domain_name = Keys::cleanName(Keys::cleanPListName(plist));
  if (!plist.isParameter("file name base")) plist.set<std::string>("file name base", "checkpoint");

  Utils::IOFileFormat_kind file_format =
    Utils::createIOFileFormat(plist, Utils::IOFileFormat_kind::HDF5);
  switch (file_format) {
  case Utils::IOFileFormat_kind::HDF5:
    return std::make_unique<OutputHDF5>(plist, comm);
  default: {
    Errors::Message msg;
    msg << "OutputFactory::createForCheckpoint: invalid \"file format\" \""
        << Utils::to_string(file_format) << "\" for writing a checkpoint file";
    throw(msg);
  }
  }
}


} // namespace OutputFactory
} // namespace Amanzi
