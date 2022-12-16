/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, others
*/

//! Utility functions for working with various file formats and names.
#ifndef AMANZI_MESH_FILE_TYPE_HH_
#define AMANZI_MESH_FILE_TYPE_HH_

#include "AmanziTypes.hh"

namespace Amanzi {
namespace AmanziMesh {

// Identifers for those file formats understood
enum struct FileFormat {
  UNKNOWN = 0, /* It's a mystery */
  EXODUS_II,   /* Exodus II format */
  NEMESIS,     /* Exodus II format partitioned by Nemesis */
  MOAB_HDF5    /* HDF5 format used by MOAB */
};

// Get the name of a particular file format
std::string
fileFormatName(const FileFormat f);

// Determine a file format from a file name, and perform sanity checks on that
// filename.
FileFormat
fileFormatFromFilename(const Comm_type& comm_, std::string name);

} // namespace AmanziMesh
} // namespace Amanzi

#endif
