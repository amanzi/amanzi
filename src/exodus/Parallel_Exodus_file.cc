// -------------------------------------------------------------
/**
 * @file   Parallel_Exodus_file.cc
 * @author William A. Perkins
 * @date Mon May  2 13:05:09 2011
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 15, 2010 by William A. Perkins
// Last Change: Mon May  2 13:05:09 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <algorithm>
#include <exception>
#include <boost/format.hpp>
#include <exodusII.h>

#include "Parallel_Exodus_file.hh"
#include "Exodus_readers.hh"


namespace Amanzi {
namespace Exodus {

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

  // create the format string that pads me with the correct number of zeros
  int ndigits = (int)floor(log10(np)) + 1;
  std::stringstream frmt;
  frmt << ".%d.%0" << ndigits << "d";

  s += boost::str(boost::format( frmt.str().c_str() ) % np % me);
  //s += boost::str(boost::format(".%d.%d") % np % me);

  std::cerr << "Process " << me << " of " << np << ": trying file " << s << std::endl;

  int ierr(0);
  ExodusError ex_all;

  try {
    my_file.reset(new Exodus_file(s.c_str()));
  } catch (const ExodusError& e) {
    std::string msg;
    msg = boost::str(boost::format("Process %d: error opening \"%s\": %s") %
                     me % s % e.what());
    std::cerr << msg << std::endl;
    ierr += 1;
    ex_all.add_data(msg.c_str());
  }

  int gerr(0);
  my_comm->SumAll(&ierr, &gerr, 1);

  if (gerr > 0) {
    Exceptions::amanzi_throw( ex_all );
  }
}

Parallel_Exodus_file::~Parallel_Exodus_file(void)
{
  // nothing to do
}

// -------------------------------------------------------------
// Parallel_Exodus_file::read_mesh
// -------------------------------------------------------------
Teuchos::RCP<AmanziMesh::Data::Data> 
Parallel_Exodus_file::read_mesh(void)
{

  my_mesh.reset(read_exodus_file(*my_file));

  // Even though element blocks are defined in all local files, if the
  // local block is empty, it will not have a cell type specified, so
  // we need get that from the other processes

  my_comm->Barrier();
  const int np(my_comm->NumProc());
  const int me(my_comm->MyPID());

  std::vector<int> byproc(np);
  int nblk(my_mesh->element_blocks());

  // check the number of blocks; should be the same on all processes

  int ierr(0);
  my_comm->GatherAll(&nblk, &byproc[0], 1);
  for (int p = 1; p < np; p++) {
    if (byproc[p] != byproc[p-1]) ierr++;
  }
  int aerr(0);
  my_comm->SumAll(&ierr, &aerr, 1);

  if (aerr) {
    std::string msg(my_basename);
    msg += ": mismatched numbers of element blocks";
    Exceptions::amanzi_throw( ExodusError(msg.c_str()) );
  }

  for (int b = 0; b < nblk; b++) {
    int mytype(my_mesh->element_block(b).element_type());
    std::vector<int> alltype(np, AmanziMesh::UNKNOWN);
    my_comm->GatherAll(&mytype, &alltype[0], 1);

    std::vector<int>::iterator junk;
    junk = std::remove(alltype.begin(), alltype.end(),
                       AmanziMesh::UNKNOWN);
    alltype.erase(junk, alltype.end());
    junk = std::unique(alltype.begin(), alltype.end());
    alltype.erase(junk, alltype.end());
    
    if (alltype.empty()) {

      // this means that all processes reported the AmanziMesh::Data::UNKNOWN
      // type for this element block; this is OK as long as it's empty
      // on all processes.

      if (my_mesh->element_block(b).num_elements() > 0) {
        std::string msg = 
          boost::str(boost::format("Process %d: %s: element block %d: block element type unknown") %
                     me % my_basename %  my_mesh->element_block(b).id());
        std::cerr << msg << std::endl;
        ierr++;
      }
    } else if (alltype.size() > 1) {
      
      // this means that at least two processes reported different,
      // not unknown, element types for this block; this is an error
      // and needs to be reported.

      std::string msg = 
        boost::str(boost::format("Process %d: %s: element block %d: block element type mismatch") %
                   me % my_basename %  my_mesh->element_block(b).id());
      std::cerr << msg << std::endl;
      ierr++;

    } else {

      AmanziMesh::Cell_type
        thetype(static_cast<AmanziMesh::Cell_type>(alltype.front()));

      if (my_mesh->element_block(b).element_type() == AmanziMesh::UNKNOWN) {
        my_mesh->element_block(b).element_type(thetype);
      }
    }
  }

  my_comm->SumAll(&ierr, &aerr, 1);
  if (aerr) {
    std::string msg = 
      boost::str(boost::format("%s: element block type errors") % my_basename);
    Exceptions::amanzi_throw( ExodusError(msg.c_str()) );
  }

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
  if (ret_val < 0) {
    std::string msg = 
      boost::str(boost::format("%s: error: cannot read element number map (%d)") %
                 my_file->filename % ret_val);
    Exceptions::amanzi_throw( ExodusError (msg.c_str()) );
  }

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
  if (ret_val < 0) {
    std::string msg = 
      boost::str(boost::format("%s: error: cannot read vertex number map (%d)") %
                 my_file->filename % ret_val);
    Exceptions::amanzi_throw( ExodusError (msg.c_str()) );
  }

  my_comm->Barrier();

  Teuchos::RCP<Epetra_Map> vmap(new Epetra_Map(-1, myvert, &gids[0], 1, *my_comm));
  
  return vmap;
}

} // namespace Exodus
} // namespace Amanzi


