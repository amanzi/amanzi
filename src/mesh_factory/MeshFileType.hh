// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: MeshFileType.hh
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 11, 2011 by William A. Perkins
// Last Change: Mon Mar 14 08:58:36 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _MeshFileType_hh_
#define _MeshFileType_hh_

#include <Epetra_Comm.h>

namespace Mesh {

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
  extern Format file_format(const Epetra_Comm& comm, const char *name);

  /// Determine, if possible, the format of the specified file
  inline Format file_format(const Epetra_Comm& comm, const std::string& name)
  {
    return file_format(comm, name.c_str());
  }

} // namespace Mesh


#endif
