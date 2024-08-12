/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! NetCDFReader: simple reader for serial reads of NetCDF files.
#include "NetCDFReader.hh"

namespace Amanzi {

NetCDFReader::NetCDFReader(const std::string& filename) : filename_(filename), file_(-1)
{
  int ierr = nc_open(filename.c_str(), NC_NOWRITE, &file_);
  if (ierr) {
    Errors::Message msg;
    msg << "NetCDFReader: error opening file \"" << filename << "\" with NC_NOWRITE access.";
    Exceptions::amanzi_throw(msg);
  }
  AMANZI_ASSERT(file_);
}


NetCDFReader::~NetCDFReader()
{
  if (file_) nc_close(file_);
}


std::pair<int,int>
NetCDFReader::findVarOrGroup_(std::string lvarname) const
{
  int ncid = file_;
  int varid = -1;
  if (lvarname[0] == '/') lvarname = lvarname.substr(1, std::string::npos);

  bool done = false;
  while (!done) {
    KeyPair name = Keys::split(lvarname, '/');
    int ncid2 = 0;
    if (name.first != "") {
      // not last entry, recurse
      int ierr = nc_inq_ncid(ncid, name.first.c_str(), &ncid2);
      if (ierr) {
        return { -1, -1 };
      } else {
        ncid = ncid2;
      }

    } else {
      // last entry... could be group or variable
      // is it a group?
      int ierr = nc_inq_ncid(ncid, name.second.c_str(), &ncid2);
      if (ierr) {
        // is it a variable?
        ierr = nc_inq_varid(ncid, name.second.c_str(), &varid);
      } else {
        ncid = ncid2;
      }
      if (ierr) return { -1, -1 };
      done = true;
    }
    lvarname = name.second;
  }
  return { ncid, varid };
}


bool
NetCDFReader::hasVariableOrGroup(const std::string& varname) const
{
  auto ncid = findVarOrGroup_(varname);
  return ncid.first >= 0 || ncid.second >= 0;
}

} // namespace Amanzi
