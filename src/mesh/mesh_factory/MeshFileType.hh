/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

#ifndef AMANZI_MESH_FILE_TYPE_HH_
#define AMANZI_MESH_FILE_TYPE_HH_

#include "AmanziTypes.hh"
#include "Teuchos_Comm.hpp"

namespace Amanzi {
namespace AmanziMesh {

  /// Identifers for those file formats understood
  enum Format { 
    UnknownFormat = 0,          /**< It's a mystery */
    ExodusII,                   /**< Exodus II format */
    Nemesis,                    /**< Exodus II format partitioned by Nemesis */
    MOABHDF5                    /**< HDF5 format used by MOAB */
  };

  /// Get the name of a particular file format
  extern std::string file_format_name(const Format& f);

  /// Determine, if possible, the format of the specified file
  extern Format file_format(Comm_ptr_type comm, const char *name);

  /// Determine, if possible, the format of the specified file
  inline Format file_format(Comm_ptr_type comm, const std::string& name)
  {
    return file_format(comm, name.c_str());
  }

} // namespace AmanziMesh
} // namespace Amanzi

#endif
