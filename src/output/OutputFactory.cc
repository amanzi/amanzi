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

  std::string output_type = plist.get<std::string>("file format", "XDMF");
  if (output_type == "XDMF" || output_type == "xdmf" || output_type == "Xdmf") {
    return std::make_unique<OutputXDMF>(plist, mesh);

  } else if (output_type == "SILO" || output_type == "silo" || output_type == "Silo") {
#ifdef HAVE_SILO
    return std::make_unique<OutputSilo>(plist, mesh);
#else
    Errors::Message msg(
      "OutputFactory: requested output type Silo was not enabled at configure time.");
    throw(msg);
#endif
  } else {
    Errors::Message msg;
    msg << "OutputFactory: Unknown output type \"" << output_type << "\"";
    throw(msg);
  }
}

std::unique_ptr<Output>
createForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm)
{
  // set the filename base based on domain if needed
  std::string domain_name = Keys::cleanName(Keys::cleanPListName(plist));
  if (!plist.isParameter("file name base")) plist.set<std::string>("file name base", "checkpoint");

  std::string output_type = plist.get<std::string>("file format", "HDF5");
  if (output_type == "HDF5") {
    return std::make_unique<OutputHDF5>(plist, comm);
  } else {
    Errors::Message msg;
    msg << "OutputFactory: Unknown output type \"" << output_type << "\"";
    throw(msg);
  }
}


} // namespace OutputFactory
} // namespace Amanzi
