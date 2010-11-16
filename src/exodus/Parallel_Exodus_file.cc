// -------------------------------------------------------------
/**
 * @file   Parallel_Exodus_file.cc
 * @author William A. Perkins
 * @date Tue Nov 16 07:13:47 2010
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 15, 2010 by William A. Perkins
// Last Change: Tue Nov 16 07:13:47 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <exception>
#include <boost/format.hpp>
#include <exodusII.h>

#include "Parallel_Exodus_file.hh"
#include "Exodus_readers.hh"


namespace ExodusII
{

// -------------------------------------------------------------
//  class Parallel_Exodus_file
// -------------------------------------------------------------

// -------------------------------------------------------------
// Parallel_Exodus_file:: constructors / destructor
// -------------------------------------------------------------
/** 
 * 
 * The files are name as @c nem_spread names them, namely @e
 * basename.N.n, where @e N is the total number of processors and @e n
 * is the local processor rank.
 * 
 * @param comm parallel environment
 * @param basename ExodusII file set base name
 */
Parallel_Exodus_file::Parallel_Exodus_file(const Epetra_Comm& comm, 
                                           const std::string& basename)
  : my_comm(comm.Clone()), my_basename(basename)
{
  const int np(my_comm->NumProc());
  const int me(my_comm->MyPID());

  std::string s(my_basename);

  s += boost::str(boost::format(".%d.%d") % np % me);

  std::cerr << "Process " << me << " of " << np << ": trying file " << s << std::endl;

  int ierr(0);
  try {
    my_file.reset(new Exodus_file(s.c_str()));
  } catch (const ExodusError& e) {
    std::string msg;
    msg = boost::str(boost::format("Process %d: error opening \"%s\": %s") %
                     me % s % e.what());
    std::cerr << msg << std::endl;
    ierr += 1;
  }

  int gerr(0);
  my_comm->SumAll(&ierr, &gerr, 1);

  if (gerr > 0) {
    throw std::runtime_error("Parallel_Exodus_file construction error");
  }
}

Parallel_Exodus_file::~Parallel_Exodus_file(void)
{
  // nothing to do
}

// -------------------------------------------------------------
// Parallel_Exodus_file::read_mesh
// -------------------------------------------------------------
Teuchos::RCP<Mesh_data::Data> 
Parallel_Exodus_file::read_mesh(void)
{

  my_mesh.reset(ExodusII::read_exodus_file(*my_file));

  return my_mesh;
}

// -------------------------------------------------------------
// Parallel_Exodus_file::cellmap
// -------------------------------------------------------------
Teuchos::RCP<Epetra_Map>
Parallel_Exodus_file::cellmap(void)
{

  if (my_mesh.is_null()) {
    read_mesh();
  }
  
  int myelem(my_mesh->parameters().num_elements_);

  std::vector<int> gids(myelem);

  int ret_val = 
    ex_get_elem_num_map(my_file->id, &gids[0]);
  if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

  my_comm->Barrier();

  Teuchos::RCP<Epetra_Map> cmap(new Epetra_Map(-1, myelem, &gids[0], 1, *my_comm));
  
  return cmap;
}

// -------------------------------------------------------------
// Parallel_Exodus_file::vertexmap
// -------------------------------------------------------------
Teuchos::RCP<Epetra_Map>
Parallel_Exodus_file::vertexmap(void)
{

  if (my_mesh.is_null()) {
    read_mesh();
  }
  
  int myvert(my_mesh->parameters().num_nodes_);

  std::vector<int> gids(myvert);

  int ret_val = 
    ex_get_node_num_map(my_file->id, &gids[0]);
  if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

  my_comm->Barrier();

  Teuchos::RCP<Epetra_Map> vmap(new Epetra_Map(-1, myvert, &gids[0], 1, *my_comm));
  
  return vmap;
}

} // close namespace ExodusII


