/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
// file: MeshFileType.cc
// -------------------------------------------------------------
/**
 * @file   MeshFileType.cc
 * @author William A. Perkins
 * @date Thu Jul 28 14:16:00 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created March 11, 2011 by William A. Perkins
// Last Change: Thu Jul 28 14:16:00 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <boost/regex.hpp>

#define BOOST_FILESYSTEM_VERSION 2
#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

#include <boost/format.hpp>


#include <fstream>

#include "MeshFileType.hh"
#include "MeshException.hh"

namespace Amanzi {
namespace AmanziMesh {

  static const std::string HDF5magic = "\211HDF\r\n\032\n";
  static const std::string NetCDFmagic1 = "CDF\001";
  static const std::string NetCDFmagic2 = "CDF\002";
  static const int magiclen = 16;

  static const boost::regex ExodusExt(".*\\.exo$", boost::regex::basic);
  static const boost::regex HDF5Ext(".*\\.h5m$", boost::regex::basic);
  static const boost::regex NemesisExt(".*\\.par$", boost::regex::basic);

  // -------------------------------------------------------------
  // file_format_name
  // -------------------------------------------------------------
  std::string
  file_format_name(const Format& f)
  {
    std::string result;
    switch (f) {
    case UnknownFormat:
    default:
      result = "Unknown";
      break;
    case ExodusII:
      result = "ExodusII";
      break;
    case MOABHDF5:
      result = "HDF5 (MOAB)";
      break;
    }
    return result;
  }

  // -------------------------------------------------------------
  // file_format
  // -------------------------------------------------------------
  /** 
   * Collective
   *
   * This routine identifies the format of the mesh file specified by @c name.  
   *
   * Currently, this is very stupid.  Format is basically only
   * determined by file name extension, as follows: 
   *
   *   - @c .exo is ::ExodusII
   *   - @c .par.N.i is ::Nemesis, where N is # cpu and i is this process id
   *   - @c .h5m is ::MOABHDF5
   * 
   * This routine also makes sure the specified file is there and is readable.
   *
   * In parallel, all processes should perform this check, even if
   * it's only one file.
   *
   * If the file exists and can be opened, this routine returns a
   * Format. If anything goes wrong, a FileMessage is thrown.
   * 
   * @param name name of file to check
   * 
   * @return format type identifier (::UnknownFormat if something goes wrong)
   */
  Format
  file_format(const Epetra_Comm& comm, const char *name) 
  {
    const int np(comm.NumProc());
    const int me(comm.MyPID());
    Format result(UnknownFormat);

    // take a guess at the format using the file name

    std::string fname(name);
    if (boost::regex_match(fname, ExodusExt)) {
      result = ExodusII;
    } else if (boost::regex_match(fname, HDF5Ext)) {
      result = MOABHDF5;
    } else if (boost::regex_match(fname, NemesisExt)) {
      result = Nemesis;
      int ndigits = (int)floor(log10(np)) + 1;
      std::string fmt = 
        boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
      fname = boost::str(boost::format(fmt) % 
                         name % comm.NumProc() % comm.MyPID());
    } else {
      result = UnknownFormat;
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
    case (ExodusII):
    case (Nemesis):
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
    case (MOABHDF5):
      fmagic.assign(buffer, HDF5magic.size());
      if (fmagic == HDF5magic) {
        ok = true;
      } else {
        e.add_data(": bad magic number, expected HDF5");
      }
      break;
    }
    if (!ok) {
      amanzi_throw(e);
    }

    return result;
  }

} // close namespace AmanziMesh
} // close namespace Amanzi
