/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

#ifndef PARALLEL_EXODUS_FILE_HH_
#define PARALLEL_EXODUS_FILE_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziMap.hh"

#include "AmanziTypes.hh"
#include "Data.hh"
#include "Exodus_file.hh"


namespace Amanzi {
namespace Exodus {

// -------------------------------------------------------------
//  class Parallel_Exodus_file
// -------------------------------------------------------------
// A reader of partitioned ExodusII files.  
/**
 * This class reads a partitioned ExodusII file set and creates a
 * Mesh_data::Data instance representing the mesh.
 *
 * The files are name as @c nem_spread names them, namely @e
 * basename.N.n, where @e N is the total number of processors and @e n
 * is the local processor rank.
 * 
 */
class Parallel_Exodus_file {
 public:

  // Default constructor.
  Parallel_Exodus_file(const Comm_ptr_type& comm, const std::string& basename);

  // Destructor
  ~Parallel_Exodus_file(void) = default;

  // Read the (local) mesh from the file
  Teuchos::RCP<AmanziMesh::Data::Data> read_mesh(void);

  // Construct a cell map for the file set (collective)
  Map_ptr_type cellmap(void);

  // Construct a vertex map for the file set (collective)
  Map_ptr_type vertexmap(void);

  Comm_ptr_type Comm() const { return comm_; }

 protected:
  Comm_ptr_type comm_;
  std::string basename_;      /**< The Exodus files' base name */

  Teuchos::RCP<Exodus_file> file_; /**< The local Exodus file description */

  Teuchos::RCP<AmanziMesh::Data::Data> mesh_; /**< The local mesh */

  Parallel_Exodus_file(const Parallel_Exodus_file& old) = delete;

  
};


} // namespace Exodus
} // namespace Amanzi

#endif
