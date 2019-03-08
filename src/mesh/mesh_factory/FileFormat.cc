/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, others
*/

//! Utility functions for working with various file formats and names.

#include <regex>
#include <fstream>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/format.hpp"
namespace bfs = boost::filesystem;

#include "FileFormat.hh"
#include "MeshException.hh"

namespace Amanzi {
namespace AmanziMesh {

static const std::string HDF5magic = "\211HDF\r\n\032\n";
static const std::string NetCDFmagic1 = "CDF\001";
static const std::string NetCDFmagic2 = "CDF\002";
static const int magiclen = 16;

static const std::regex ExodusExt(".*\\.exo$");
static const std::regex HDF5Ext(".*\\.h5m$");
static const std::regex NemesisExt(".*\\.par$");

// -------------------------------------------------------------
// file_format_name
// -------------------------------------------------------------
std::string fileFormatName(const FileFormat f)
{
  std::string result;
  switch (f) {
  case FileFormat::UNKNOWN:
  default:
    result = "Unknown";
    break;
  case FileFormat::EXODUS_II:
    result = "ExodusII";
    break;
  case FileFormat::NEMESIS:
    result = "ExodusII";
    break;
  case FileFormat::MOAB_HDF5:
    result = "HDF5 (Framework::MOAB)";
    break;
  }
  return result;
}

// -------------------------------------------------------------
// file_format
// -------------------------------------------------------------
/*
 Collective

 This routine identifies the format of the mesh file specified by @c name.  

 Currently, this is very stupid.  Format is basically only
 determined by file name extension, as follows: 

   - .exo is ::ExodusII
   - .par.N.i is ::Nemesis, where N is # cpu and i is this process id
   - .h5m is ::MOAB_HDF5
 
 This routine also makes sure the specified file is there and is readable.

 In parallel, all processes should perform this check, even if
 it's only one file.

 If the file exists and can be opened, this routine returns a
 Format. If anything goes wrong, an exception is thrown.
 */
FileFormat fileFormatFromFilename(const Comm_type& comm, std::string fname) 
{
  const int np(comm.NumProc());
  const int me(comm.MyPID());

  // take a guess at the format using the file name
  auto result = FileFormat::UNKNOWN;

  if (std::regex_match(fname, ExodusExt)) {
    result = FileFormat::EXODUS_II;
  } else if (std::regex_match(fname, HDF5Ext)) {
    result = FileFormat::MOAB_HDF5;
  } else if (std::regex_match(fname, NemesisExt)) {
    result = FileFormat::NEMESIS;
    int ndigits = (int)floor(log10(np)) + 1;
    std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
    fname = boost::str(boost::format(fmt) % 
                       fname % np % me);
  }

  // check to see if there is actually a file first
  FileMessage e(fname.c_str());
  bfs::path p(fname);
  if (!bfs::exists(p)) {
    e.add_data(": path not found");
    amanzi_throw(e);
  }
    
  // check the file's magic number
  std::ifstream s(fname.c_str(), std::ios::binary);
  if (s.fail()) {
    e.add_data(": cannot open");
    amanzi_throw(e);
  }
    
  char buffer[magiclen];
  s.read(buffer, magiclen);
  s.close();

  std::string fmagic;
  bool ok(false);
  switch (result) {
  case (FileFormat::UNKNOWN):
      Exceptions::amanzi_throw(e);
      break;
  case (FileFormat::EXODUS_II):
  case (FileFormat::NEMESIS):
    fmagic.assign(buffer, NetCDFmagic1.size());
    if (fmagic == NetCDFmagic1) { 
      ok = true;
    } 
    fmagic.assign(buffer, NetCDFmagic2.size());
    if (fmagic == NetCDFmagic2) { 
      ok = true;
    } 
    if (!ok) {
      e.add_data(": bad magic number, expected NetCDF");
      fmagic.assign(buffer, HDF5magic.size());
      if (fmagic == HDF5magic) {
        ok = true;
      } else {
        e.add_data(", expected HDF5");
      }
    } 
    break;
  case (FileFormat::MOAB_HDF5):
    fmagic.assign(buffer, HDF5magic.size());
    if (fmagic == HDF5magic) {
      ok = true;
    } else {
      e.add_data(": bad magic number, expected HDF5");
    }
    break;
  }
  if (!ok) {
    Exceptions::amanzi_throw(e);
  }

  return result;
}

} // namespace AmanziMesh
} // namespace Amanzi
