// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   Parallel_Exodus_file.hh
 * @author William A. Perkins
 * @date Tue Nov 16 10:09:55 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 15, 2010 by William A. Perkins
// Last Change: Tue Nov 16 10:09:55 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _Parallel_Exodus_file_hh_
#define _Parallel_Exodus_file_hh_

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "Data.hh"
#include "Exodus_file.hh"


namespace ExodusII
{

// -------------------------------------------------------------
//  class Parallel_Exodus_file
// -------------------------------------------------------------
/// A reader of partitioned ExodusII files.  
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
protected:

  
  Teuchos::RCP<Epetra_Comm> my_comm; /**< The parallel environment */
  std::string my_basename;      /**< The Exodus files' base name */

  Teuchos::RCP<Exodus_file> my_file; /**< The local Exodus file description */

  Teuchos::RCP<Mesh_data::Data> my_mesh; /**< The local mesh */

  /// Protected, unimplemented, copy constructor to avoid unwanted copies.
  Parallel_Exodus_file(const Parallel_Exodus_file& old);

public:

  /// Default constructor.
  Parallel_Exodus_file(const Epetra_Comm& comm, const std::string& basename);

  /// Destructor
  ~Parallel_Exodus_file(void);

  /// Get the parallel environment
  const Epetra_Comm& comm() { return *my_comm; }

  /// Read the (local) mesh from the file
  Teuchos::RCP<Mesh_data::Data> read_mesh(void);

  /// Construct a cell map for the file set (collective)
  Teuchos::RCP<Epetra_Map> cellmap(void);

  /// Construct a vertex map for the file set (collective)
  Teuchos::RCP<Epetra_Map> vertexmap(void);
};


} // close namespace ExodusII

#endif
